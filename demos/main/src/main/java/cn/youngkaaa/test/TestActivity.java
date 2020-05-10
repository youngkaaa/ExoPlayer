package cn.youngkaaa.test;

import android.net.Uri;
import android.os.Environment;
import android.os.Handler;
import android.os.Looper;
import androidx.appcompat.app.AppCompatActivity;
import android.os.Bundle;
import com.google.android.exoplayer2.SimpleExoPlayer;
import com.google.android.exoplayer2.demo.R;
import com.google.android.exoplayer2.source.MediaSource;
import com.google.android.exoplayer2.source.ProgressiveMediaSource;
import com.google.android.exoplayer2.ui.PlayerView;
import com.google.android.exoplayer2.upstream.DefaultDataSourceFactory;
import java.io.File;

public class TestActivity extends AppCompatActivity {

  @Override
  protected void onCreate(Bundle savedInstanceState) {
    super.onCreate(savedInstanceState);
    setContentView(R.layout.activity_test);


    PlayerView playerView = findViewById(R.id.player_view);

    SimpleExoPlayer player = new SimpleExoPlayer.Builder(this)
        .build();

    String url = "https://fanyiapp.cdn.bcebos.com/ugc/video/captain_marvel.mp4";
    Uri uri = Uri.parse(url);

    DefaultDataSourceFactory dataSourceFactory = new DefaultDataSourceFactory(this, "test");

    MediaSource mediaSource = new ProgressiveMediaSource.Factory(dataSourceFactory)
        .createMediaSource(uri);

    playerView.setPlayer(player);
    player.setPlayWhenReady(true);
    player.prepare(mediaSource);
  }
}
