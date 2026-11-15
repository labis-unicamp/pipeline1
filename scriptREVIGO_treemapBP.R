# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.
revigo.names <- c("term_ID","description","frequency","uniqueness","dispensability","representative");
revigo.data <- rbind(
c("GO:0005201","extracellular matrix structural constituent",0.184,0.918,0.0,"extracellular matrix structural constituent"),
c("GO:0005515","protein binding",13.416,0.777,0.0,"protein binding"),
c("GO:0005576","extracellular region",2.961,0.939,0.0,"extracellular region"),
c("GO:0005615","extracellular space",0.967,0.920,0.0,"extracellular space"),
c("GO:0005796","Golgi lumen",0.057,0.930,0.0,"Golgi lumen"),
c("GO:0005886","plasma membrane",7.096,0.828,0.0,"plasma membrane"),
c("GO:0031012","extracellular matrix",0.595,0.892,0.0,"extracellular matrix"),
c("GO:0046872","metal ion binding",4.796,0.843,0.0,"metal ion binding"),
c("GO:0070062","extracellular exosome",1.700,0.863,0.0,"extracellular exosome"),
c("GO:0070701","mucus layer",0.007,0.939,0.0,"mucus layer"),
c("GO:0000226","microtubule cytoskeleton organization",0.354,0.812,0.0,"microtubule cytoskeleton organization"),
c("GO:0001725","stress fiber",0.051,0.875,0.0,"stress fiber"),
c("GO:0003779","actin binding",0.671,0.818,0.0,"actin binding"),
c("GO:0005178","integrin binding",0.113,0.892,0.0,"integrin binding"),
c("GO:0005198","structural molecule activity",0.896,0.918,0.0,"structural molecule activity"),
c("GO:0005509","calcium ion binding",1.032,0.834,0.0,"calcium ion binding"),
c("GO:0005604","basement membrane",0.124,0.917,0.0,"basement membrane"),
c("GO:0005634","nucleus",6.618,0.870,0.0,"nucleus"),
c("GO:0005635","nuclear envelope",0.229,0.881,0.0,"nuclear envelope"),
c("GO:0005654","nucleoplasm",1.049,0.865,0.0,"nucleoplasm"),
c("GO:0005737","cytoplasm",7.515,0.883,0.0,"cytoplasm"),
c("GO:0005783","endoplasmic reticulum",1.178,0.861,0.0,"endoplasmic reticulum"),
c("GO:0005789","endoplasmic reticulum membrane",0.552,0.842,0.0,"endoplasmic reticulum membrane"),
c("GO:0005829","cytosol",3.184,0.885,0.0,"cytosol"),
c("GO:0005856","cytoskeleton",1.234,0.809,0.0,"cytoskeleton"),
c("GO:0005874","microtubule",0.478,0.843,0.0,"microtubule"),
c("GO:0005882","intermediate filament",0.131,0.875,0.0,"intermediate filament"),
c("GO:0005925","focal adhesion",0.267,0.847,0.0,"focal adhesion"),
c("GO:0005938","cell cortex",0.099,0.894,0.0,"cell cortex"),
c("GO:0007010","cytoskeleton organization",1.367,0.787,0.0,"cytoskeleton organization"),
c("GO:0007155","cell adhesion",1.332,0.840,0.0,"cell adhesion"),
c("GO:0007229","integrin-mediated signaling pathway",0.156,0.849,0.0,"integrin-mediated signaling pathway"),
c("GO:0008017","microtubule binding",0.324,0.847,0.0,"microtubule binding"),
c("GO:0008090","retrograde axonal transport",0.010,0.897,0.0,"retrograde axonal transport"),
c("GO:0009611","response to wounding",0.577,0.918,0.0,"response to wounding"),
c("GO:0015629","actin cytoskeleton",0.464,0.838,0.0,"actin cytoskeleton"),
c("GO:0015630","microtubule cytoskeleton",0.656,0.828,0.0,"microtubule cytoskeleton"),
c("GO:0016020","membrane",10.136,0.918,0.0,"membrane"),
c("GO:0030011","maintenance of cell polarity",0.029,0.870,0.0,"maintenance of cell polarity"),
c("GO:0030018","Z disc",0.049,0.886,0.0,"Z disc"),
c("GO:0030056","hemidesmosome",0.024,0.895,0.0,"hemidesmosome"),
c("GO:0030424","axon",0.331,0.847,0.0,"axon"),
c("GO:0031252","cell leading edge",0.194,0.851,0.0,"cell leading edge"),
c("GO:0031410","cytoplasmic vesicle",1.258,0.897,0.0,"cytoplasmic vesicle"),
c("GO:0031673","H zone",0.002,0.927,0.0,"H zone"),
c("GO:0035371","microtubule plus-end",0.006,0.896,0.0,"microtubule plus-end"),
c("GO:0042060","wound healing",0.248,0.880,0.0,"wound healing"),
c("GO:0042803","protein homodimerization activity",0.652,0.866,0.0,"protein homodimerization activity"),
c("GO:0042995","cell projection",1.213,0.847,0.0,"cell projection"),
c("GO:0045104","intermediate filament cytoskeleton organization",0.056,0.855,0.0,"intermediate filament cytoskeleton organization"),
c("GO:0045111","intermediate filament cytoskeleton",0.146,0.855,0.0,"intermediate filament cytoskeleton"),
c("GO:0048870","cell motility",0.652,0.847,0.0,"cell motility"),
c("GO:0051010","microtubule plus-end binding",0.022,0.870,0.0,"microtubule plus-end binding"),
c("GO:0070161","anchoring junction",0.324,0.861,0.0,"anchoring junction"),
c("GO:0071944","cell periphery",4.854,0.879,0.0,"cell periphery"),
c("GO:0099513","polymeric cytoskeletal fiber",0.860,0.830,0.0,"polymeric cytoskeletal fiber"),
c("GO:1904115","axon cytoplasm",0.026,0.895,0.0,"axon cytoplasm"),
c("GO:0005085","guanyl-nucleotide exchange factor activity",0.337,0.885,0.0,"guanyl-nucleotide exchange factor activity"),
c("GO:0005096","GTPase activator activity",0.275,0.867,0.0,"GTPase activator activity"),
c("GO:0005768","endosome",0.444,0.890,0.0,"endosome"),
c("GO:0006897","endocytosis",0.509,0.841,0.0,"endocytosis"),
c("GO:0016192","vesicle-mediated transport",1.126,0.831,0.0,"vesicle-mediated transport"),
c("GO:0030139","endocytic vesicle",0.119,0.888,0.0,"endocytic vesicle"),
c("GO:0031267","small GTPase binding",0.558,0.824,0.0,"small GTPase binding"),
c("GO:0043087","regulation of GTPase activity",0.631,0.810,0.0,"regulation of GTPase activity"),
c("GO:0045296","cadherin binding",0.173,0.865,0.0,"cadherin binding"),
c("GO:0051223","regulation of protein transport",0.423,0.815,0.0,"regulation of protein transport")
);

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "uniqueness",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()
