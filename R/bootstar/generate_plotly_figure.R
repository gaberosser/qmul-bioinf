library(plotly)

trace1 <- list(
  x = c(
    273.72523786414,
    320.728035056893,
    303.782633873609,
    285.083117289354,
    280.630303975773,
    277.152537748482,
    297.998457828263,
    241.150045780902,
    251.537188370809,
    255.113767354,
    268.158379349657,
    262.531307250501
  ),
  
  y = c(
    206.541665872345,
    85.7923110460854,
    157.218874018819,
    192.113714488529,
    56.7192257421416,
    27.874625039053,
    54.6509285786413,
    233.811378608119,
    6.44783153443546,
    17.0570415160178,
    95.6360953212528,
    127.643209511883
  ),
  
  z = c(
    -97.5357101225432,
    -10.5446378763575,
    -90.5058608939255,
    -93.5725969319984,
    -110.741584111114,
    -67.7933138360432,
    -40.5878433591888,
    -82.7538609755062,
    -26.2922179642214,
    -20.9351272461744,
    -113.722533589883,
    -111.616046891126
  ),
  
  error_x = list(color = "rgba(102,194,165,1)"),
  
  error_y = list(color = "rgba(102,194,165,1)"),
  
  #line = list(color = "rgba(102,194,165,1)"),
  
  marker = list(
    color = "rgba(102,194,165,1)",
    
    line = list(color = "rgba(102,194,165,1)"),
    
    symbol = "square"
    
  ),
  
  mode = "markers",
  
  name = "Astrocyte<br />Boot et al.",
  
  textfont = list(color = "rgba(102,194,165,1)"),
  
  type = "scatter3d",
  
  xsrc = "JBOOT:4:586b2a",
  
  ysrc = "JBOOT:4:f67534",
  
  zsrc = "JBOOT:4:bbcb13"
  
)

trace2 <- list(
  x = c(247.066929579377),
  
  y = c(-156.150539802935),
  
  z = c(50.8304473329749),
  
  error_x = list(color = "rgba(102,194,165,1)"),
  
  error_y = list(color = "rgba(102,194,165,1)"),
  
  #line = list(color = "rgba(102,194,165,1)"),
  
  marker = list(
    color = "rgba(102,194,165,1)",
    
    line = list(color = "rgba(102,194,165,1)"),
    
    symbol = "diamond"
    
  ),
  
  mode = "markers",
  
  name = "Astrocyte<br />Zhou et al.",
  
  textfont = list(color = "rgba(102,194,165,1)"),
  
  type = "scatter3d",
  
  xsrc = "JBOOT:4:4a43e3",
  
  ysrc = "JBOOT:4:83d0b8",
  
  zsrc = "JBOOT:4:e890de"
  
)

trace3 <- list(
  x = c(
    316.243103340574,
    254.580460378129,
    264.415810271412,
    285.345393502805
  ),
  
  y = c(
    -290.484726049253,
    -271.938272408372,
    -266.844643814666,
    -279.158252004525
  ),
  
  z = c(
    184.421538115362,
    163.123432858298,
    159.579246630945,
    175.075580835047
  ),
  
  error_x = list(color = "rgba(252,141,98,1)"),
  
  error_y = list(color = "rgba(252,141,98,1)"),
  
  #line = list(color = "rgba(252,141,98,1)"),
  
  marker = list(
    color = "rgba(252,141,98,1)",
    
    line = list(color = "rgba(252,141,98,1)"),
    
    symbol = "square"
    
  ),
  
  mode = "markers",
  
  name = "Fibroblast<br />Boot et al.",
  
  textfont = list(color = "rgba(252,141,98,1)"),
  
  type = "scatter3d",
  
  xsrc = "JBOOT:4:baeaf3",
  
  ysrc = "JBOOT:4:d6a582",
  
  zsrc = "JBOOT:4:9075d0"
  
)

trace4 <- list(
  x = c(
    -186.821433580896,
    -201.52656890801,
    -151.104655334583,
    -186.821433580896,
    -201.52656890801
  ),
  
  y = c(
    -254.986080453448,
    -232.751226144491,
    -253.727858778677,
    -254.986080453447,
    -232.751226144491
  ),
  
  z = c(
    -117.601253529034,
    -132.968254176372,
    -123.316281818476,
    -117.601253529033,
    -132.968254176372
  ),
  
  error_x = list(color = "rgba(141,160,203,1)"),
  
  error_y = list(color = "rgba(141,160,203,1)"),
  
  #line = list(color = "rgba(141,160,203,1)"),
  
  marker = list(
    color = "rgba(141,160,203,1)",
    
    line = list(color = "rgba(141,160,203,1)"),
    
    symbol = "square"
    
  ),
  
  mode = "markers",
  
  name = "iPSC<br />Boot et al.",
  
  textfont = list(color = "rgba(141,160,203,1)"),
  
  type = "scatter3d",
  
  xsrc = "JBOOT:4:b883b5",
  
  ysrc = "JBOOT:4:9e56ea",
  
  zsrc = "JBOOT:4:d90264"
  
)

trace5 <- list(
  x = c(
    -103.191592927696,
    -95.6242480219896,
    -109.978351615294,
    -199.001350048069,
    -179.796979250497,
    -206.820925362793,
    -169.084523504348,
    -200.702413517775,
    -177.940488817313,
    -171.832227556774,
    -123.955831916214,
    -115.90554772611
  ),
  
  y = c(
    -2.44223654493369,
    -53.596579392515,
    -38.8915973954671,
    -19.959898783023,
    -8.74913379450649,
    -48.9932067440702,
    91.4474450122622,
    -31.4476916780501,
    12.9594562709851,
    -67.7297475998256,
    -23.7693236254687,
    -10.7294328561878
  ),
  
  z = c(
    -44.0630184276683,
    -32.2755017696872,
    -24.3007176909636,
    -54.7931676832564,
    -60.3002878786001,
    -55.1222560002955,
    -22.0894800286934,
    -59.0875890931794,
    -58.5806866568478,
    -68.4657464773382,
    -45.4702636522639,
    -60.3503524896042
  ),
  
  error_x = list(color = "rgba(231,138,195,1)"),
  
  error_y = list(color = "rgba(231,138,195,1)"),
  
  #line = list(color = "rgba(231,138,195,1)"),
  
  marker = list(
    color = "rgba(231,138,195,1)",
    
    line = list(color = "rgba(231,138,195,1)"),
    
    symbol = "square"
    
  ),
  
  mode = "markers",
  
  name = "NSC<br />Boot et al.",
  
  textfont = list(color = "rgba(231,138,195,1)"),
  
  type = "scatter3d",
  
  xsrc = "JBOOT:4:84dd9a",
  
  ysrc = "JBOOT:4:324eb5",
  
  zsrc = "JBOOT:4:1a4a4a"
  
)

trace6 <- list(
  x = c(-165.524858907456),
  
  y = c(174.398050972073),
  
  z = c(49.3965258691586),
  
  error_x = list(color = "rgba(231,138,195,1)"),
  
  error_y = list(color = "rgba(231,138,195,1)"),
  
  #line = list(color = "rgba(231,138,195,1)"),
  
  marker = list(
    color = "rgba(231,138,195,1)",
    
    line = list(color = "rgba(231,138,195,1)"),
    
    symbol = "cross"
    
  ),
  
  mode = "markers",
  
  name = "NSC<br />Gibco",
  
  textfont = list(color = "rgba(231,138,195,1)"),
  
  type = "scatter3d",
  
  xsrc = "JBOOT:4:5fc3d8",
  
  ysrc = "JBOOT:4:aa105a",
  
  zsrc = "JBOOT:4:e0bc71"
  
)

trace7 <- list(
  x = c(
    -204.898221951228,
    -155.78784468263,
    -183.266918557762,
    -165.381308563845,
    -5.86174249901382,
    -160.433979982332,
    -153.715341040543,
    -170.947862473688,
    -185.238469573689,
    -193.893663476594,
    -159.85853560434,
    -55.7978306357778
  ),
  
  y = c(
    170.70109525695,
    177.4865755531,
    107.839526680377,
    121.924867093218,
    82.754642583814,
    111.158485729078,
    92.1956841919525,
    139.689919366881,
    93.5234769597434,
    150.566206670484,
    71.8070180655242,
    53.7307968515442
  ),
  
  z = c(
    103.308985625292,
    112.4523838401,
    90.4321291087757,
    118.71616653846,
    147.85040421082,
    83.4195096937308,
    75.2491711868059,
    91.7456818745547,
    83.6612453053046,
    83.190181691888,
    97.9331509745846,
    146.619946761651
  ),
  
  error_x = list(color = "rgba(166,216,84,1)"),
  
  error_y = list(color = "rgba(166,216,84,1)"),
  
  #line = list(color = "rgba(166,216,84,1)"),
  
  marker = list(
    color = "rgba(166,216,84,1)",
    
    line = list(color = "rgba(166,216,84,1)"),
    
    symbol = "square"
    
  ),
  
  mode = "markers",
  
  name = "OPC<br />Boot et al.",
  
  textfont = list(color = "rgba(166,216,84,1)"),
  
  type = "scatter3d",
  
  xsrc = "JBOOT:4:c8587f",
  
  ysrc = "JBOOT:4:fb82a9",
  
  zsrc = "JBOOT:4:3c0b2c"
  
)

trace8 <- list(
  x = c(56.999009711484),
  
  y = c(-113.602394066953),
  
  z = c(58.9499704220125),
  
  error_x = list(color = "rgba(166,216,84,1)"),
  
  error_y = list(color = "rgba(166,216,84,1)"),
  
  #line = list(color = "rgba(166,216,84,1)"),
  
  marker = list(
    color = "rgba(166,216,84,1)",
    
    line = list(color = "rgba(166,216,84,1)"),
    
    symbol = "diamond"
    
  ),
  
  mode = "markers",
  
  name = "OPC<br />Zhou et al.",
  
  textfont = list(color = "rgba(166,216,84,1)"),
  
  type = "scatter3d",
  
  xsrc = "JBOOT:4:1b71e4",
  
  ysrc = "JBOOT:4:5f4118",
  
  zsrc = "JBOOT:4:d715b4"
  
)

data <-
  list(trace1, trace2, trace3, trace4, trace5, trace6, trace7, trace8)

layout <- list(
  hovermode = "closest",
  
  margin = list(
    r = 10,
    
    t = 25,
    
    b = 40,
    
    l = 60
    
  ),
  
  scene = list(
    xaxis = list(title = "PCA1"),
    
    yaxis = list(title = "PCA2"),
    
    zaxis = list(title = "PCA3")
    
  ),
  
  showlegend = TRUE
  
)

p <- plot_ly()

p <-
  add_trace(
    p,
    x = trace1$x,
    y = trace1$y,
    z = trace1$z,
    error_x = trace1$error_x,
    error_y = trace1$error_y,
    line = trace1$line,
    marker = trace1$marker,
    mode = trace1$mode,
    name = trace1$name,
    textfont = trace1$textfont,
    type = trace1$type,
    xsrc = trace1$xsrc,
    ysrc = trace1$ysrc,
    zsrc = trace1$zsrc
  )

p <-
  add_trace(
    p,
    x = trace2$x,
    y = trace2$y,
    z = trace2$z,
    error_x = trace2$error_x,
    error_y = trace2$error_y,
    line = trace2$line,
    marker = trace2$marker,
    mode = trace2$mode,
    name = trace2$name,
    textfont = trace2$textfont,
    type = trace2$type,
    xsrc = trace2$xsrc,
    ysrc = trace2$ysrc,
    zsrc = trace2$zsrc
  )

p <-
  add_trace(
    p,
    x = trace3$x,
    y = trace3$y,
    z = trace3$z,
    error_x = trace3$error_x,
    error_y = trace3$error_y,
    line = trace3$line,
    marker = trace3$marker,
    mode = trace3$mode,
    name = trace3$name,
    textfont = trace3$textfont,
    type = trace3$type,
    xsrc = trace3$xsrc,
    ysrc = trace3$ysrc,
    zsrc = trace3$zsrc
  )

p <-
  add_trace(
    p,
    x = trace4$x,
    y = trace4$y,
    z = trace4$z,
    error_x = trace4$error_x,
    error_y = trace4$error_y,
    line = trace4$line,
    marker = trace4$marker,
    mode = trace4$mode,
    name = trace4$name,
    textfont = trace4$textfont,
    type = trace4$type,
    xsrc = trace4$xsrc,
    ysrc = trace4$ysrc,
    zsrc = trace4$zsrc
  )

p <-
  add_trace(
    p,
    x = trace5$x,
    y = trace5$y,
    z = trace5$z,
    error_x = trace5$error_x,
    error_y = trace5$error_y,
    line = trace5$line,
    marker = trace5$marker,
    mode = trace5$mode,
    name = trace5$name,
    textfont = trace5$textfont,
    type = trace5$type,
    xsrc = trace5$xsrc,
    ysrc = trace5$ysrc,
    zsrc = trace5$zsrc
  )

p <-
  add_trace(
    p,
    x = trace6$x,
    y = trace6$y,
    z = trace6$z,
    error_x = trace6$error_x,
    error_y = trace6$error_y,
    line = trace6$line,
    marker = trace6$marker,
    mode = trace6$mode,
    name = trace6$name,
    textfont = trace6$textfont,
    type = trace6$type,
    xsrc = trace6$xsrc,
    ysrc = trace6$ysrc,
    zsrc = trace6$zsrc
  )

p <-
  add_trace(
    p,
    x = trace7$x,
    y = trace7$y,
    z = trace7$z,
    error_x = trace7$error_x,
    error_y = trace7$error_y,
    line = trace7$line,
    marker = trace7$marker,
    mode = trace7$mode,
    name = trace7$name,
    textfont = trace7$textfont,
    type = trace7$type,
    xsrc = trace7$xsrc,
    ysrc = trace7$ysrc,
    zsrc = trace7$zsrc
  )

p <-
  add_trace(
    p,
    x = trace8$x,
    y = trace8$y,
    z = trace8$z,
    error_x = trace8$error_x,
    error_y = trace8$error_y,
    line = trace8$line,
    marker = trace8$marker,
    mode = trace8$mode,
    name = trace8$name,
    textfont = trace8$textfont,
    type = trace8$type,
    xsrc = trace8$xsrc,
    ysrc = trace8$ysrc,
    zsrc = trace8$zsrc
  )

p <-
  layout(
    p,
    hovermode = layout$hovermode,
    margin = layout$margin,
    scene = layout$scene,
    showlegend = layout$showlegend
  )



packageVersion('plotly')


library(RSelenium)
rD <- RSelenium::rsDriver(browser = 'firefox')
export(
  p = p,
  file = 'my_figure.svg',
  selenium = rD
)
