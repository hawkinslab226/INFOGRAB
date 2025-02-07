##colors  vector
mycolors1 <- c(
  `Kidney`="#43009A",
  `Trachea`="#990099",
  `BCell`="#0DFAFA",
  `Spleen_T_Cells`="#13B7B7",
  `Bursa`="#004C99",
  `Thymus`="#2685E4",
  `d0_macrophage`="#F02B6D",
  `d3_macrophage`="#FF0091",
  `d6_macrophage`="#F572BC",
  `Macrophage.lung`="#D06AAA",
  `Monocyte.blood,`="#91155B",
  `Ileum`="#98D55A",
  `Jejunum`="#4C9900",
  `Proximal.Cecum`="#CCFF99",
  `Dark.meat`="#A1122A",
  `White.meat`="#FFCCCC",
  `Isthmus.od2`="#FF8000",
  `Magnum.od1`="#DE5100",
  `Shell.Gland.od3`="#FFC78E",
  `Ovary`="#CCCC00"
)

#function colors
mycols1 <- function(...) {
  cols <- c(...)
  if (is.null(cols))
    return (mycolors1)
  mycolors1[cols]
}

###To request the color tissue use the names
color_names= names(mycolors1)
color_names