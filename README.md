This repository contains some simple [R
Shiny](https://shiny.rstudio.com/) applications that I developed to
illustrate topics I cover in my graduate course on [Population
Genetics](http://darwin.eeb.uconn.edu/uncommon-ground/eeb348/). Each
of the applications is hosted at shinyapps.io. You'll find them
embedded in various lecture details on my course website. Below you'll
find a list of direct links to make access a bit more convenient.

* [Estimating allele frequencies with the EM algorithm](https://keholsinger.shinyapps.io/EM-algorithm-for-allele-frequencies/)
* [Genetics of geographically structured populations](https://keholsinger.shinyapps.io/F-statistics/)
* [Viability selection at one locus with two alleles](https://keholsinger.shinyapps.io/Viability-selection/)
* [Allele frequency changes with genetic drift](https://keholsinger.shinyapps.io/Genetic-Drift/)
* [Allele frequency changes with genetic drift and mutation](https://keholsinger.shinyapps.io/Drift-mutation/)
* [The genealogy of the coalescent at one locus in one population](https://keholsinger.shinyapps.io/coalescent/)
* [The genealogy of the coalescent at one locus in two populations with migration and mutation](https://github.com/kholsinger/PopGen-Shiny/blob/master/Coalescent-structured-population/app.R) Note: Unlike the links above, this is a link to the R Shiny source. I use <tt>ggtree()</tt> to display and color the coalescent tree. Unfortunately, I can't get <tt>shinyapps.io</tt> to install the application right now. You'll need to download the source to your hard drive (Click "Raw" and save it someplace convenient), boot up <tt>R</tt>, load Shiny (<tt>library(shiny)</tt>), make sure your working directory is the same as the directory where you saved <tt>app.R</tt>, and <tt>runApp()</tt>. I hope to get this fixed, but no promises on when that will happen.