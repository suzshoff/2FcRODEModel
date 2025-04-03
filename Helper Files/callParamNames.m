function paramNames = paramNames(fcr1, fcr2)

paramNames = strings(1,31);

paramNames(1) = "k-IgG1-Antigen-On";
paramNames(2) = "k-IgG1-Antigen-Off";
paramNames(3) = "k-IgG2-Antigen-On";
paramNames(4) = "k-IgG2-Antigen-Off";
paramNames(5) = "k-IgG3-Antigen-On";
paramNames(6) = "k-IgG3-Antigen-Off";
paramNames(7) = "k-IgG4-Antigen-On";
paramNames(8) = "k-IgG4-Antigen-Off";

paramNames(9) = strcat("k-IgG1-", fcr1, "-On");
paramNames(10) = strcat("k-IgG1-", fcr1, "-Off");
paramNames(11) = strcat("k-IgG2-", fcr1, "-On");
paramNames(12) = strcat("k-IgG2-", fcr1, "-Off");
paramNames(13) = strcat("k-IgG3-", fcr1, "-On");
paramNames(14) = strcat("k-IgG3-", fcr1, "-Off");
paramNames(15) = strcat("k-IgG4-", fcr1, "-On");
paramNames(16) = strcat("k-IgG4-", fcr1, "-Off");

paramNames(17) = strcat("k-IgG1-", fcr2, "-On");
paramNames(18) = strcat("k-IgG1-", fcr2, "-Off");
paramNames(19) = strcat("k-IgG2-", fcr2, "-On");
paramNames(20) = strcat("k-IgG2-", fcr2, "-Off");
paramNames(21) = strcat("k-IgG3-", fcr2, "-On");
paramNames(22) = strcat("k-IgG3-", fcr2, "-Off");
paramNames(23) = strcat("k-IgG4-", fcr2, "-On");
paramNames(24) = strcat("k-IgG4-", fcr2, "-Off");

paramNames(25) = "IgG1";
paramNames(26) = "IgG2";
paramNames(27) = "IgG3";
paramNames(28) = "IgG4";
paramNames(29) = "Antigen";
paramNames(30) = strcat("FcR", fcr1);
paramNames(31) = strcat("FcR", fcr2);
