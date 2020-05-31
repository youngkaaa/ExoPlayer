/*
 * Copyright (C) 2016 The Android Open Source Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package com.google.android.exoplayer2.extractor.mp4;

import androidx.annotation.IntDef;
import com.google.android.exoplayer2.C;
import com.google.android.exoplayer2.Format;
import com.google.android.exoplayer2.ParserException;
import com.google.android.exoplayer2.audio.Ac4Util;
import com.google.android.exoplayer2.extractor.Extractor;
import com.google.android.exoplayer2.extractor.ExtractorInput;
import com.google.android.exoplayer2.extractor.ExtractorOutput;
import com.google.android.exoplayer2.extractor.ExtractorsFactory;
import com.google.android.exoplayer2.extractor.GaplessInfoHolder;
import com.google.android.exoplayer2.extractor.PositionHolder;
import com.google.android.exoplayer2.extractor.SeekMap;
import com.google.android.exoplayer2.extractor.SeekPoint;
import com.google.android.exoplayer2.extractor.TrackOutput;
import com.google.android.exoplayer2.extractor.mp4.Atom.ContainerAtom;
import com.google.android.exoplayer2.metadata.Metadata;
import com.google.android.exoplayer2.util.Assertions;
import com.google.android.exoplayer2.util.MimeTypes;
import com.google.android.exoplayer2.util.NalUnitUtil;
import com.google.android.exoplayer2.util.ParsableByteArray;
import java.io.IOException;
import java.lang.annotation.Documented;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;

/**
 * Extracts data from the MP4 container format.
 */
public final class Mp4Extractor implements Extractor, SeekMap {

  /**
   * Factory for {@link Mp4Extractor} instances.
   */
  public static final ExtractorsFactory FACTORY = () -> new Extractor[]{new Mp4Extractor()};

  /**
   * Flags controlling the behavior of the extractor. Possible flag value is {@link
   * #FLAG_WORKAROUND_IGNORE_EDIT_LISTS}.
   */
  @Documented
  @Retention(RetentionPolicy.SOURCE)
  @IntDef(
      flag = true,
      value = {FLAG_WORKAROUND_IGNORE_EDIT_LISTS})
  public @interface Flags {

  }

  /**
   * Flag to ignore any edit lists in the stream.
   */
  public static final int FLAG_WORKAROUND_IGNORE_EDIT_LISTS = 1;

  /**
   * Parser states.
   */
  @Documented
  @Retention(RetentionPolicy.SOURCE)
  @IntDef({STATE_READING_ATOM_HEADER, STATE_READING_ATOM_PAYLOAD, STATE_READING_SAMPLE})
  private @interface State {

  }

  private static final int STATE_READING_ATOM_HEADER = 0;
  private static final int STATE_READING_ATOM_PAYLOAD = 1;
  private static final int STATE_READING_SAMPLE = 2;

  /**
   * Brand stored in the ftyp atom for QuickTime media.
   */
  private static final int BRAND_QUICKTIME = 0x71742020; // ftyp: qt

  /**
   * When seeking within the source, if the offset is greater than or equal to this value (or the
   * offset is negative), the source will be reloaded.
   */
  private static final long RELOAD_MINIMUM_SEEK_DISTANCE = 256 * 1024;

  /**
   * For poorly interleaved streams, the maximum byte difference one track is allowed to be read
   * ahead before the source will be reloaded at a new position to read another track.
   */
  private static final long MAXIMUM_READ_AHEAD_BYTES_STREAM = 10 * 1024 * 1024;

  private final @Flags
  int flags;

  // Temporary arrays.
  private final ParsableByteArray nalStartCode;
  private final ParsableByteArray nalLength;
  private final ParsableByteArray scratch;

  private final ParsableByteArray atomHeader;
  private final ArrayDeque<ContainerAtom> containerAtoms;

  @State
  private int parserState;
  private int atomType;
  private long atomSize;
  private int atomHeaderBytesRead; // 当前已读取的box(atom)header的长度 每次enterReadingAtomHeaderState也就是重新读取下一个box(atom)时都会置为0
  private ParsableByteArray atomData;

  private int sampleTrackIndex;
  private int sampleBytesRead;
  private int sampleBytesWritten;
  private int sampleCurrentNalBytesRemaining;

  // Extractor outputs.
  private ExtractorOutput extractorOutput;
  private Mp4Track[] tracks;
  private long[][] accumulatedSampleSizes;
  private int firstVideoTrackIndex; // 第一个视频轨道索引值
  private long durationUs; // 所有轨道中最大的时长，也就是该视频文件的时长
  private boolean isQuickTime; // 当前mp4 ftyp中brand是否是quickTime的

  /**
   * Creates a new extractor for unfragmented MP4 streams.
   */
  public Mp4Extractor() {
    this(0);
  }

  /**
   * Creates a new extractor for unfragmented MP4 streams, using the specified flags to control the
   * extractor's behavior.
   *
   * @param flags Flags that control the extractor's behavior.
   */
  public Mp4Extractor(@Flags int flags) {
    this.flags = flags;
    atomHeader = new ParsableByteArray(Atom.LONG_HEADER_SIZE);
    containerAtoms = new ArrayDeque<>();
    nalStartCode = new ParsableByteArray(NalUnitUtil.NAL_START_CODE);
    nalLength = new ParsableByteArray(4);
    scratch = new ParsableByteArray();
    sampleTrackIndex = C.INDEX_UNSET;
  }

  @Override
  public boolean sniff(ExtractorInput input) throws IOException, InterruptedException {
    return Sniffer.sniffUnfragmented(input);
  }

  @Override
  public void init(ExtractorOutput output) {
    extractorOutput = output;
  }

  @Override
  public void seek(long position, long timeUs) {
    containerAtoms.clear();
    atomHeaderBytesRead = 0;
    sampleTrackIndex = C.INDEX_UNSET;
    sampleBytesRead = 0;
    sampleBytesWritten = 0;
    sampleCurrentNalBytesRemaining = 0;
    if (position == 0) {
      enterReadingAtomHeaderState();
    } else if (tracks != null) {
      updateSampleIndices(timeUs);
    }
  }

  @Override
  public void release() {
    // Do nothing
  }

  @Override
  public int read(ExtractorInput input, PositionHolder seekPosition)
      throws IOException, InterruptedException {
    while (true) {
      switch (parserState) {
        // 读取box的header，也就是固定的8字节长度 size(4)+type(4)
        case STATE_READING_ATOM_HEADER:
          if (!readAtomHeader(input)) {
            return RESULT_END_OF_INPUT;
          }
          break;
        // 读取header之后的body
        case STATE_READING_ATOM_PAYLOAD:
          if (readAtomPayload(input, seekPosition)) {
            // 需要外部seek一下，seek的位置在positionHolder.position中保存
            return RESULT_SEEK;
          }
          break;
        case STATE_READING_SAMPLE:
          return readSample(input, seekPosition);
        default:
          throw new IllegalStateException();
      }
    }
  }

  // SeekMap implementation.

  @Override
  public boolean isSeekable() {
    return true;
  }

  @Override
  public long getDurationUs() {
    return durationUs;
  }

  @Override
  public SeekPoints getSeekPoints(long timeUs) {
    if (tracks.length == 0) {
      return new SeekPoints(SeekPoint.START);
    }

    long firstTimeUs;
    long firstOffset;
    long secondTimeUs = C.TIME_UNSET;
    long secondOffset = C.POSITION_UNSET;

    // If we have a video track, use it to establish one or two seek points.
    if (firstVideoTrackIndex != C.INDEX_UNSET) {
      TrackSampleTable sampleTable = tracks[firstVideoTrackIndex].sampleTable;
      int sampleIndex = getSynchronizationSampleIndex(sampleTable, timeUs);
      if (sampleIndex == C.INDEX_UNSET) {
        return new SeekPoints(SeekPoint.START);
      }
      long sampleTimeUs = sampleTable.timestampsUs[sampleIndex];
      firstTimeUs = sampleTimeUs;
      firstOffset = sampleTable.offsets[sampleIndex];
      if (sampleTimeUs < timeUs && sampleIndex < sampleTable.sampleCount - 1) {
        int secondSampleIndex = sampleTable.getIndexOfLaterOrEqualSynchronizationSample(timeUs);
        if (secondSampleIndex != C.INDEX_UNSET && secondSampleIndex != sampleIndex) {
          secondTimeUs = sampleTable.timestampsUs[secondSampleIndex];
          secondOffset = sampleTable.offsets[secondSampleIndex];
        }
      }
    } else {
      firstTimeUs = timeUs;
      firstOffset = Long.MAX_VALUE;
    }

    // Take into account other tracks.
    for (int i = 0; i < tracks.length; i++) {
      if (i != firstVideoTrackIndex) {
        TrackSampleTable sampleTable = tracks[i].sampleTable;
        firstOffset = maybeAdjustSeekOffset(sampleTable, firstTimeUs, firstOffset);
        if (secondTimeUs != C.TIME_UNSET) {
          secondOffset = maybeAdjustSeekOffset(sampleTable, secondTimeUs, secondOffset);
        }
      }
    }

    SeekPoint firstSeekPoint = new SeekPoint(firstTimeUs, firstOffset);
    if (secondTimeUs == C.TIME_UNSET) {
      return new SeekPoints(firstSeekPoint);
    } else {
      SeekPoint secondSeekPoint = new SeekPoint(secondTimeUs, secondOffset);
      return new SeekPoints(firstSeekPoint, secondSeekPoint);
    }
  }

  // Private methods.

  /**
   * 进入到读取Header的模式
   */
  private void enterReadingAtomHeaderState() {
    parserState = STATE_READING_ATOM_HEADER;
    atomHeaderBytesRead = 0;
  }

  /**
   * 读取Mp4中每个box(atom)中的前8个字节。包括:length(4)+type(4)
   */
  private boolean readAtomHeader(ExtractorInput input) throws IOException, InterruptedException {
    // 每次切换进入到STATE_READING_ATOM_HEADER时，都会将其置为0.表示开始读取新的box了
    if (atomHeaderBytesRead == 0) {
      // Read the standard length atom header.
      // 读取Atom.HEADER_SIZE(8个)字节到atomHeader.data中去。读取失败则返回false。读取成功则往下走
      if (!input.readFully(atomHeader.data, 0, Atom.HEADER_SIZE, true)) {
        return false;
      }

      atomHeaderBytesRead = Atom.HEADER_SIZE; // 读取成功的话，那么header应该是读取了Atom.HEADER_SIZE也就是8个字节
      atomHeader.setPosition(0); // atomHeader位置重置，因为里面data更新了，要从data的头部开始读取
      atomSize = atomHeader.readUnsignedInt(); // 当前box长度大小 字节数，也就是4个字节，atomHeader里面position往后挪动4
      atomType = atomHeader.readInt(); // 当前box的类型，也就是4个字节，atomHeader里面position往后挪动4
    }
    // size=1时表示其大小是large size。真正的size得在largesize域上得到。（实际上只有“mdat”类型的box才有可能用到large size。）
    if (atomSize == Atom.DEFINES_LARGE_SIZE) {
      // Read the large size.
      int headerBytesRemaining = Atom.LONG_HEADER_SIZE - Atom.HEADER_SIZE;
      input.readFully(atomHeader.data, Atom.HEADER_SIZE, headerBytesRemaining);
      atomHeaderBytesRead += headerBytesRemaining;
      atomSize = atomHeader.readUnsignedLongToLong();
    } else if (atomSize == Atom.EXTENDS_TO_END_SIZE) {
      // 如果size为0，表示该box为文件的最后一个box，文件结尾即为该box结尾。（同样只存在于“mdat”类型的box中。）
      // The atom extends to the end of the file. Note that if the atom is within a container we can
      // work out its size even if the input length is unknown.
      long endPosition = input.getLength();
      if (endPosition == C.LENGTH_UNSET && !containerAtoms.isEmpty()) {
        endPosition = containerAtoms.peek().endPosition;
      }
      if (endPosition != C.LENGTH_UNSET) {
        atomSize = endPosition - input.getPosition() + atomHeaderBytesRead;
      }
    }
    // atomSize不合法。肯定要大于等于8，等于8的话表示没有后面的body data
    if (atomSize < atomHeaderBytesRead) {
      throw new ParserException("Atom size less than header length (unsupported).");
    }
    // 判断该box类型是否是container box类型
    if (shouldParseContainerAtom(atomType)) {
      // 本box container的结束position
      long endPosition =
          input.getPosition() + atomSize - atomHeaderBytesRead;
      // 针对meta做单独处理
      if (atomSize != atomHeaderBytesRead && atomType == Atom.TYPE_meta) {
        maybeSkipRemainingMetaAtomHeaderBytes(input);
      }
      // 新增加一个box container实例，添加到containerAtoms中(第一个位置)去。
      // 保存该container的endPosition，用于后续判断该box container是否读取完毕
      containerAtoms.push(new ContainerAtom(atomType, endPosition));
      // 如果该box的长度就是8也就是刚刚读取的box header，那么body没有，就直接end
      if (atomSize == atomHeaderBytesRead) {
        processAtomEnded(endPosition);
      } else {
        // Start reading the first child atom.
        // 切换为STATE_READING_ATOM_HEADER，准备读取该box container中的子container
        enterReadingAtomHeaderState();
      }
    } else if (shouldParseLeafAtom(atomType)) { // 单纯的box(atom)类型，不能包括子box(atom)
      // We don't support parsing of leaf atoms that define extended atom sizes, or that have
      // lengths greater than Integer.MAX_VALUE.
      Assertions.checkState(atomHeaderBytesRead == Atom.HEADER_SIZE);
      Assertions.checkState(atomSize <= Integer.MAX_VALUE);
      // 新创建一个atomData，保存一个box(atom)，后续在读取该box body 的readAtomPayload方法中会用到
      atomData = new ParsableByteArray((int) atomSize);
      // 把读取到的前8字节的header复制到新的atomData中去。
      System.arraycopy(atomHeader.data, 0, atomData.data, 0, Atom.HEADER_SIZE);
      // 读取该box(atom)之后剩余的body data
      parserState = STATE_READING_ATOM_PAYLOAD;
    } else { // 其他没用的box类型,比如type=free，mdat的，那么不保存，走一遍STATE_READING_ATOM_PAYLOAD把后面的body跳过
      atomData = null;
      parserState = STATE_READING_ATOM_PAYLOAD;
    }

    return true;
  }

  /**
   * Processes the atom payload. If {@link #atomData} is null and the size is at or above the
   * threshold {@link #RELOAD_MINIMUM_SEEK_DISTANCE}, {@code true} is returned and the caller should
   * restart loading at the position in {@code positionHolder}. Otherwise, the atom is
   * read/skipped.
   */
  private boolean readAtomPayload(ExtractorInput input, PositionHolder positionHolder)
      throws IOException, InterruptedException {
    // atomSize:该box的总长度；atomHeaderBytesRead：已读取的头部固定8字节长度；atomPayloadSize：该box的body大小，也就是还需要读取的大小
    long atomPayloadSize = atomSize - atomHeaderBytesRead;
    // 计算出该box相对于input的结束位置。也就是该box的结束位置
    long atomEndPosition = input.getPosition() + atomPayloadSize;
    // 是否需要seek来跳过该body。有可能该box body很长，而且该box暂时用不到，那么就需要seek来跳过。
    boolean seekRequired = false;
    // 在readAtomHeader中判断该box(atom)需要处理，那么atomData会在readAtomHeader中赋值
    if (atomData != null) {
      // 从input中读取剩余的数据，读取完成后input中的position会更新
      input.readFully(atomData.data, atomHeaderBytesRead, (int) atomPayloadSize);
      // 判读当前box(atom)类型是ftyp File Type Box，那么判断刚读取的body中内容，判断是否是quickTime类型的
      // 读取ftyp body，判断是否是quick time类型的视频，不是的话就是普通的isom
      if (atomType == Atom.TYPE_ftyp) {
        isQuickTime = processFtypAtom(atomData);
      } else if (!containerAtoms.isEmpty()) {
        // containerAtoms不为空，表示前一次在解析header时，遇到了一个box container，
        // 此时解析的是该container中的子box 所以要把该box加到该container中去
        containerAtoms.peek().add(new Atom.LeafAtom(atomType, atomData));
      }
    } else {
      // We don't need the data. Skip or seek, depending on how large the atom is.
      // 该box body我们不需要，那么判断怎么跳过，如果body长的话就seek来跳过
      // 如果要跳过的数据小于RELOAD_MINIMUM_SEEK_DISTANCE，就通过skipFully来跳过
      if (atomPayloadSize < RELOAD_MINIMUM_SEEK_DISTANCE) {
        input.skipFully((int) atomPayloadSize);
      } else { // 需要跳过的长度太大了，就把seekRequired置为true，用重新请求的方式来跳过
        // 更新最终seek的position，这个等到后面seek时会用到。
        positionHolder.position = input.getPosition() + atomPayloadSize;
        seekRequired = true;
      }
    }
    // 处理完一个box(atom)就走一遍这个方法。把当前该box的结束位置传入进去。
    processAtomEnded(atomEndPosition);
    return seekRequired && parserState != STATE_READING_SAMPLE;
  }

  /**
   * 每读完一个box就调用执行本方法，用于判断moov box container是否已读完
   *
   * @param atomEndPosition 当前读完的box的结束位置
   */
  private void processAtomEnded(long atomEndPosition) throws ParserException {
    // 如果containerAtoms不为空，标明现在是在读取box container中的子box；
    // 然后判断当前读完的box的结束位置和之前存储的container box的endPosition是否一样，一样的话表示该box container已经处理完毕。
    while (!containerAtoms.isEmpty() && containerAtoms.peek().endPosition == atomEndPosition) {
      // 从containerAtoms中拿出该box container
      Atom.ContainerAtom containerAtom = containerAtoms.pop();
      // 如果该box container是moov box container，那么标明moov box container 读取完毕了
      // 此时调用processMoovAtom方法来处理解析moov box container，并且切到STATE_READING_SAMPLE模式准备读取采样数据了
      if (containerAtom.type == Atom.TYPE_moov) {
        // We've reached the end of the moov atom. Process it and prepare to read samples.
        processMoovAtom(containerAtom); // 处理moov box。处理完成后轨道就有了，并且外部onPrepare会被回调
        containerAtoms.clear();
        parserState = STATE_READING_SAMPLE; // 进入读取采样的模式了。
      } else if (!containerAtoms.isEmpty()) {
        // 该box不是moov，但是处理完毕了，此时如果它外部还有box container，说明它是该box container的子box，所以将其加入进去
        containerAtoms.peek().add(containerAtom);
      }
    }
    // 如果moov读取完了，会将其切到STATE_READING_SAMPLE模式；
    // 如果读取完的不是moov，那么此时肯定不是STATE_READING_SAMPLE，所以开始读取下一个box，所以进入到STATE_READING_ATOM_HEADER模式
    if (parserState != STATE_READING_SAMPLE) {
      enterReadingAtomHeaderState();
    }
  }

  /**
   * Updates the stored track metadata to reflect the contents of the specified moov atom. moov box
   * container读取完毕了。此时开始处理该box container
   */
  private void processMoovAtom(ContainerAtom moov) throws ParserException {
    int firstVideoTrackIndex = C.INDEX_UNSET; // 第一个视频轨道的索引
    long durationUs = C.TIME_UNSET;
    List<Mp4Track> tracks = new ArrayList<>();

    // Process metadata.
    Metadata udtaMetadata = null;
    GaplessInfoHolder gaplessInfoHolder = new GaplessInfoHolder();
    Atom.LeafAtom udta = moov.getLeafAtomOfType(Atom.TYPE_udta); // udta:用户数据
    if (udta != null) {
      udtaMetadata = AtomParsers.parseUdta(udta, isQuickTime);
      if (udtaMetadata != null) {
        gaplessInfoHolder.setFromMetadata(udtaMetadata);
      }
    }
    Metadata mdtaMetadata = null;
    Atom.ContainerAtom meta = moov.getContainerAtomOfType(Atom.TYPE_meta); // meta
    if (meta != null) {
      mdtaMetadata = AtomParsers.parseMdtaFromMeta(meta);
    }

    boolean ignoreEditLists = (flags & FLAG_WORKAROUND_IGNORE_EDIT_LISTS) != 0;
    // 解析每个track中的信息。
    ArrayList<TrackSampleTable> trackSampleTables =
        getTrackSampleTables(moov, gaplessInfoHolder, ignoreEditLists);

    int trackCount = trackSampleTables.size(); // 轨道数
    for (int i = 0; i < trackCount; i++) {
      TrackSampleTable trackSampleTable = trackSampleTables.get(i);
      Track track = trackSampleTable.track;
      long trackDurationUs =
          track.durationUs != C.TIME_UNSET ? track.durationUs : trackSampleTable.durationUs;
      // 取轨道中最大的时长作为最终的时长
      durationUs = Math.max(durationUs, trackDurationUs);
      Mp4Track mp4Track = new Mp4Track(track, trackSampleTable,
          extractorOutput.track(i, track.type));

      // Each sample has up to three bytes of overhead for the start code that replaces its length.
      // Allow ten source samples per output sample, like the platform extractor.
      int maxInputSize = trackSampleTable.maximumSize + 3 * 10;
      Format format = track.format.copyWithMaxInputSize(maxInputSize);
      if (track.type == C.TRACK_TYPE_VIDEO
          && trackDurationUs > 0
          && trackSampleTable.sampleCount > 1) {
        // 计算帧率。视频采样数/视频时长(秒)
        float frameRate = trackSampleTable.sampleCount / (trackDurationUs / 1000000f);
        format = format.copyWithFrameRate(frameRate); // 更新帧率
      }
      format =
          MetadataUtil.getFormatWithMetadata(
              track.type, format, udtaMetadata, mdtaMetadata, gaplessInfoHolder);
      mp4Track.trackOutput.format(format);

      if (track.type == C.TRACK_TYPE_VIDEO && firstVideoTrackIndex == C.INDEX_UNSET) {
        firstVideoTrackIndex = tracks.size(); // 保存第一个视频轨道的索引
      }
      tracks.add(mp4Track);
    }
    this.firstVideoTrackIndex = firstVideoTrackIndex;
    this.durationUs = durationUs;
    this.tracks = tracks.toArray(new Mp4Track[0]);
    accumulatedSampleSizes = calculateAccumulatedSampleSizes(this.tracks);

    extractorOutput.endTracks(); // track解析完毕
    extractorOutput.seekMap(this);
  }

  /**
   * 解析出每个轨道。以及每个轨道中的采样数据表
   *
   * @param moov
   * @param gaplessInfoHolder
   * @param ignoreEditLists
   * @return
   * @throws ParserException
   */
  private ArrayList<TrackSampleTable> getTrackSampleTables(
      ContainerAtom moov, GaplessInfoHolder gaplessInfoHolder, boolean ignoreEditLists)
      throws ParserException {
    ArrayList<TrackSampleTable> trackSampleTables = new ArrayList<>();
    for (int i = 0; i < moov.containerChildren.size(); i++) {
      Atom.ContainerAtom atom = moov.containerChildren.get(i);
      // 不是trak box的就跳过
      if (atom.type != Atom.TYPE_trak) {
        continue;
      }
      // 解析该track
      Track track =
          AtomParsers.parseTrak(
              atom,
              moov.getLeafAtomOfType(Atom.TYPE_mvhd),
              /* duration= */ C.TIME_UNSET,
              /* drmInitData= */ null,
              ignoreEditLists,
              isQuickTime);
      if (track == null) {
        continue;
      }
      Atom.ContainerAtom stblAtom =
          atom.getContainerAtomOfType(Atom.TYPE_mdia)
              .getContainerAtomOfType(Atom.TYPE_minf)
              .getContainerAtomOfType(Atom.TYPE_stbl);
      TrackSampleTable trackSampleTable = AtomParsers.parseStbl(track, stblAtom, gaplessInfoHolder);
      if (trackSampleTable.sampleCount == 0) {
        continue;
      }
      trackSampleTables.add(trackSampleTable); // 该track解析完毕，加入到列表中
    }
    return trackSampleTables;
  }

  /**
   * Attempts to extract the next sample in the current mdat atom for the specified track.
   * <p>
   * Returns {@link #RESULT_SEEK} if the source should be reloaded from the position in {@code
   * positionHolder}.
   * <p>
   * Returns {@link #RESULT_END_OF_INPUT} if no samples are left. Otherwise, returns {@link
   * #RESULT_CONTINUE}.
   *
   * @param input          The {@link ExtractorInput} from which to read data.
   * @param positionHolder If {@link #RESULT_SEEK} is returned, this holder is updated to hold the
   *                       position of the required data.
   * @return One of the {@code RESULT_*} flags in {@link Extractor}.
   * @throws IOException          If an error occurs reading from the input.
   * @throws InterruptedException If the thread is interrupted.
   */
  private int readSample(ExtractorInput input, PositionHolder positionHolder)
      throws IOException, InterruptedException {
    long inputPosition = input.getPosition();
    if (sampleTrackIndex == C.INDEX_UNSET) {
      sampleTrackIndex = getTrackIndexOfNextReadSample(inputPosition);
      if (sampleTrackIndex == C.INDEX_UNSET) {
        return RESULT_END_OF_INPUT;
      }
    }
    Mp4Track track = tracks[sampleTrackIndex];
    TrackOutput trackOutput = track.trackOutput;
    int sampleIndex = track.sampleIndex;
    long position = track.sampleTable.offsets[sampleIndex];
    int sampleSize = track.sampleTable.sizes[sampleIndex];
    long skipAmount = position - inputPosition + sampleBytesRead;
    if (skipAmount < 0 || skipAmount >= RELOAD_MINIMUM_SEEK_DISTANCE) {
      positionHolder.position = position;
      return RESULT_SEEK;
    }
    if (track.track.sampleTransformation == Track.TRANSFORMATION_CEA608_CDAT) {
      // The sample information is contained in a cdat atom. The header must be discarded for
      // committing.
      skipAmount += Atom.HEADER_SIZE;
      sampleSize -= Atom.HEADER_SIZE;
    }
    input.skipFully((int) skipAmount);
    if (track.track.nalUnitLengthFieldLength != 0) {
      // Zero the top three bytes of the array that we'll use to decode nal unit lengths, in case
      // they're only 1 or 2 bytes long.
      byte[] nalLengthData = nalLength.data;
      nalLengthData[0] = 0;
      nalLengthData[1] = 0;
      nalLengthData[2] = 0;
      int nalUnitLengthFieldLength = track.track.nalUnitLengthFieldLength;
      int nalUnitLengthFieldLengthDiff = 4 - track.track.nalUnitLengthFieldLength;
      // NAL units are length delimited, but the decoder requires start code delimited units.
      // Loop until we've written the sample to the track output, replacing length delimiters with
      // start codes as we encounter them.
      while (sampleBytesWritten < sampleSize) {
        if (sampleCurrentNalBytesRemaining == 0) {
          // Read the NAL length so that we know where we find the next one.
          input.readFully(nalLengthData, nalUnitLengthFieldLengthDiff, nalUnitLengthFieldLength);
          sampleBytesRead += nalUnitLengthFieldLength;
          nalLength.setPosition(0);
          int nalLengthInt = nalLength.readInt();
          if (nalLengthInt < 0) {
            throw new ParserException("Invalid NAL length");
          }
          sampleCurrentNalBytesRemaining = nalLengthInt;
          // Write a start code for the current NAL unit.
          nalStartCode.setPosition(0);
          trackOutput.sampleData(nalStartCode, 4);
          sampleBytesWritten += 4;
          sampleSize += nalUnitLengthFieldLengthDiff;
        } else {
          // Write the payload of the NAL unit.
          int writtenBytes = trackOutput.sampleData(input, sampleCurrentNalBytesRemaining, false);
          sampleBytesRead += writtenBytes;
          sampleBytesWritten += writtenBytes;
          sampleCurrentNalBytesRemaining -= writtenBytes;
        }
      }
    } else {
      if (MimeTypes.AUDIO_AC4.equals(track.track.format.sampleMimeType)) {
        if (sampleBytesWritten == 0) {
          Ac4Util.getAc4SampleHeader(sampleSize, scratch);
          trackOutput.sampleData(scratch, Ac4Util.SAMPLE_HEADER_SIZE);
          sampleBytesWritten += Ac4Util.SAMPLE_HEADER_SIZE;
        }
        sampleSize += Ac4Util.SAMPLE_HEADER_SIZE;
      }
      while (sampleBytesWritten < sampleSize) {
        int writtenBytes = trackOutput.sampleData(input, sampleSize - sampleBytesWritten, false);
        sampleBytesRead += writtenBytes;
        sampleBytesWritten += writtenBytes;
        sampleCurrentNalBytesRemaining -= writtenBytes;
      }
    }
    trackOutput.sampleMetadata(track.sampleTable.timestampsUs[sampleIndex],
        track.sampleTable.flags[sampleIndex], sampleSize, 0, null);
    track.sampleIndex++;
    sampleTrackIndex = C.INDEX_UNSET;
    sampleBytesRead = 0;
    sampleBytesWritten = 0;
    sampleCurrentNalBytesRemaining = 0;
    return RESULT_CONTINUE;
  }

  /**
   * Returns the index of the track that contains the next sample to be read, or {@link
   * C#INDEX_UNSET} if no samples remain.
   *
   * <p>The preferred choice is the sample with the smallest offset not requiring a source reload,
   * or if not available the sample with the smallest overall offset to avoid subsequent source
   * reloads.
   *
   * <p>To deal with poor sample interleaving, we also check whether the required memory to catch
   * up with the next logical sample (based on sample time) exceeds {@link
   * #MAXIMUM_READ_AHEAD_BYTES_STREAM}. If this is the case, we continue with this sample even
   * though it may require a source reload.
   */
  private int getTrackIndexOfNextReadSample(long inputPosition) {
    long preferredSkipAmount = Long.MAX_VALUE;
    boolean preferredRequiresReload = true;
    int preferredTrackIndex = C.INDEX_UNSET;
    long preferredAccumulatedBytes = Long.MAX_VALUE;
    long minAccumulatedBytes = Long.MAX_VALUE;
    boolean minAccumulatedBytesRequiresReload = true;
    int minAccumulatedBytesTrackIndex = C.INDEX_UNSET;
    for (int trackIndex = 0; trackIndex < tracks.length; trackIndex++) {
      Mp4Track track = tracks[trackIndex];
      int sampleIndex = track.sampleIndex;
      if (sampleIndex == track.sampleTable.sampleCount) {
        continue;
      }
      long sampleOffset = track.sampleTable.offsets[sampleIndex];
      long sampleAccumulatedBytes = accumulatedSampleSizes[trackIndex][sampleIndex];
      long skipAmount = sampleOffset - inputPosition;
      boolean requiresReload = skipAmount < 0 || skipAmount >= RELOAD_MINIMUM_SEEK_DISTANCE;
      if ((!requiresReload && preferredRequiresReload)
          || (requiresReload == preferredRequiresReload && skipAmount < preferredSkipAmount)) {
        preferredRequiresReload = requiresReload;
        preferredSkipAmount = skipAmount;
        preferredTrackIndex = trackIndex;
        preferredAccumulatedBytes = sampleAccumulatedBytes;
      }
      if (sampleAccumulatedBytes < minAccumulatedBytes) {
        minAccumulatedBytes = sampleAccumulatedBytes;
        minAccumulatedBytesRequiresReload = requiresReload;
        minAccumulatedBytesTrackIndex = trackIndex;
      }
    }
    return minAccumulatedBytes == Long.MAX_VALUE
        || !minAccumulatedBytesRequiresReload
        || preferredAccumulatedBytes < minAccumulatedBytes + MAXIMUM_READ_AHEAD_BYTES_STREAM
        ? preferredTrackIndex
        : minAccumulatedBytesTrackIndex;
  }

  /**
   * Updates every track's sample index to point its latest sync sample before/at {@code timeUs}.
   */
  private void updateSampleIndices(long timeUs) {
    for (Mp4Track track : tracks) {
      TrackSampleTable sampleTable = track.sampleTable;
      int sampleIndex = sampleTable.getIndexOfEarlierOrEqualSynchronizationSample(timeUs);
      if (sampleIndex == C.INDEX_UNSET) {
        // Handle the case where the requested time is before the first synchronization sample.
        sampleIndex = sampleTable.getIndexOfLaterOrEqualSynchronizationSample(timeUs);
      }
      track.sampleIndex = sampleIndex;
    }
  }

  /**
   * Possibly skips the version and flags fields (1+3 byte) of a full meta atom of the {@code
   * input}.
   *
   * <p>Atoms of type {@link Atom#TYPE_meta} are defined to be full atoms which have four
   * additional bytes for a version and a flags field (see 4.2 'Object Structure' in ISO/IEC
   * 14496-12:2005). QuickTime do not have such a full box structure. Since some of these files are
   * encoded wrongly, we can't rely on the file type though. Instead we must check the 8 bytes after
   * the common header bytes ourselves.
   */
  private void maybeSkipRemainingMetaAtomHeaderBytes(ExtractorInput input)
      throws IOException, InterruptedException {
    scratch.reset(8);
    // Peek the next 8 bytes which can be either
    // (iso) [1 byte version + 3 bytes flags][4 byte size of next atom]
    // (qt)  [4 byte size of next atom      ][4 byte hdlr atom type   ]
    // In case of (iso) we need to skip the next 4 bytes.
    input.peekFully(scratch.data, 0, 8);
    scratch.skipBytes(4);
    if (scratch.readInt() == Atom.TYPE_hdlr) {
      input.resetPeekPosition();
    } else {
      input.skipFully(4);
    }
  }

  /**
   * For each sample of each track, calculates accumulated size of all samples which need to be read
   * before this sample can be used.
   */
  private static long[][] calculateAccumulatedSampleSizes(Mp4Track[] tracks) {
    long[][] accumulatedSampleSizes = new long[tracks.length][];
    int[] nextSampleIndex = new int[tracks.length];
    long[] nextSampleTimesUs = new long[tracks.length];
    boolean[] tracksFinished = new boolean[tracks.length];
    for (int i = 0; i < tracks.length; i++) {
      accumulatedSampleSizes[i] = new long[tracks[i].sampleTable.sampleCount];
      nextSampleTimesUs[i] = tracks[i].sampleTable.timestampsUs[0];
    }
    long accumulatedSampleSize = 0;
    int finishedTracks = 0;
    while (finishedTracks < tracks.length) {
      long minTimeUs = Long.MAX_VALUE;
      int minTimeTrackIndex = -1;
      for (int i = 0; i < tracks.length; i++) {
        if (!tracksFinished[i] && nextSampleTimesUs[i] <= minTimeUs) {
          minTimeTrackIndex = i;
          minTimeUs = nextSampleTimesUs[i];
        }
      }
      int trackSampleIndex = nextSampleIndex[minTimeTrackIndex];
      accumulatedSampleSizes[minTimeTrackIndex][trackSampleIndex] = accumulatedSampleSize;
      accumulatedSampleSize += tracks[minTimeTrackIndex].sampleTable.sizes[trackSampleIndex];
      nextSampleIndex[minTimeTrackIndex] = ++trackSampleIndex;
      if (trackSampleIndex < accumulatedSampleSizes[minTimeTrackIndex].length) {
        nextSampleTimesUs[minTimeTrackIndex] =
            tracks[minTimeTrackIndex].sampleTable.timestampsUs[trackSampleIndex];
      } else {
        tracksFinished[minTimeTrackIndex] = true;
        finishedTracks++;
      }
    }
    return accumulatedSampleSizes;
  }

  /**
   * Adjusts a seek point offset to take into account the track with the given {@code sampleTable},
   * for a given {@code seekTimeUs}.
   *
   * @param sampleTable The sample table to use.
   * @param seekTimeUs  The seek time in microseconds.
   * @param offset      The current offset.
   * @return The adjusted offset.
   */
  private static long maybeAdjustSeekOffset(
      TrackSampleTable sampleTable, long seekTimeUs, long offset) {
    int sampleIndex = getSynchronizationSampleIndex(sampleTable, seekTimeUs);
    if (sampleIndex == C.INDEX_UNSET) {
      return offset;
    }
    long sampleOffset = sampleTable.offsets[sampleIndex];
    return Math.min(sampleOffset, offset);
  }

  /**
   * Returns the index of the synchronization sample before or at {@code timeUs}, or the index of
   * the first synchronization sample if located after {@code timeUs}, or {@link C#INDEX_UNSET} if
   * there are no synchronization samples in the table.
   *
   * @param sampleTable The sample table in which to locate a synchronization sample.
   * @param timeUs      A time in microseconds.
   * @return The index of the synchronization sample before or at {@code timeUs}, or the index of
   * the first synchronization sample if located after {@code timeUs}, or {@link C#INDEX_UNSET} if
   * there are no synchronization samples in the table.
   */
  private static int getSynchronizationSampleIndex(TrackSampleTable sampleTable, long timeUs) {
    int sampleIndex = sampleTable.getIndexOfEarlierOrEqualSynchronizationSample(timeUs);
    if (sampleIndex == C.INDEX_UNSET) {
      // Handle the case where the requested time is before the first synchronization sample.
      sampleIndex = sampleTable.getIndexOfLaterOrEqualSynchronizationSample(timeUs);
    }
    return sampleIndex;
  }

  /**
   * Process an ftyp atom to determine whether the media is QuickTime.
   *
   * @param atomData The ftyp atom data.
   * @return Whether the media is QuickTime.
   */
  private static boolean processFtypAtom(ParsableByteArray atomData) {
    atomData.setPosition(Atom.HEADER_SIZE); // 先跳过前8个固定的header
    int majorBrand = atomData.readInt(); // major brand， ftyp body中四个字节长度
    if (majorBrand == BRAND_QUICKTIME) { // 判断其是否是qt 也就是quickTime的
      return true;
    }
    atomData.skipBytes(4); // minor_version 不是的话就往后继续找
    while (atomData.bytesLeft() > 0) {
      if (atomData.readInt() == BRAND_QUICKTIME) {
        return true;
      }
    }
    return false;
  }

  /**
   * Returns whether the extractor should decode a leaf atom with type {@code atom}.
   */
  private static boolean shouldParseLeafAtom(int atom) {
    return atom == Atom.TYPE_mdhd
        || atom == Atom.TYPE_mvhd
        || atom == Atom.TYPE_hdlr
        || atom == Atom.TYPE_stsd
        || atom == Atom.TYPE_stts
        || atom == Atom.TYPE_stss
        || atom == Atom.TYPE_ctts
        || atom == Atom.TYPE_elst
        || atom == Atom.TYPE_stsc
        || atom == Atom.TYPE_stsz
        || atom == Atom.TYPE_stz2
        || atom == Atom.TYPE_stco
        || atom == Atom.TYPE_co64
        || atom == Atom.TYPE_tkhd
        || atom == Atom.TYPE_ftyp
        || atom == Atom.TYPE_udta
        || atom == Atom.TYPE_keys
        || atom == Atom.TYPE_ilst;
  }

  /**
   * Returns whether the extractor should decode a container atom with type {@code atom}.
   * 这些类型的box(atom)是container box类型，也就是其里面可以包括子box(atom)
   */
  private static boolean shouldParseContainerAtom(int atom) {
    return atom == Atom.TYPE_moov
        || atom == Atom.TYPE_trak
        || atom == Atom.TYPE_mdia
        || atom == Atom.TYPE_minf
        || atom == Atom.TYPE_stbl
        || atom == Atom.TYPE_edts
        || atom == Atom.TYPE_meta;
  }

  private static final class Mp4Track {

    public final Track track;
    public final TrackSampleTable sampleTable;
    public final TrackOutput trackOutput;

    public int sampleIndex;

    public Mp4Track(Track track, TrackSampleTable sampleTable, TrackOutput trackOutput) {
      this.track = track;
      this.sampleTable = sampleTable;
      this.trackOutput = trackOutput;
    }

  }

}
