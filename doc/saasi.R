## ----include = FALSE----------------------------------------------------------
 knitr::opts_chunk$set(
   collapse = TRUE,
   comment = "#>",
   fig.path = "man/figures/README-",
   out.width = "100%"
 )

## ----bird order data----------------------------------------------------------
library(ape)
library(saasi)
data(bird.orders)

# set some random discrete characters to the bird tree.
x <- factor(c(rep(1, 5), rep(2, 18)))

bird.orders$tip.state <- x

## ----speciation and extinction------------------------------------------------
estimates <- mle_lm(bird.orders,lambda = .2, mu = .05, psi = .1,lower = c(0.001,0.001), upper = c(0.5,0.5))

estimates

## ----tranistion---------------------------------------------------------------
ace_qij <- extract_ace_q(ace(x, bird.orders, type = "d"))

## ----saasi--------------------------------------------------------------------
# set up the parameters
pars <- data.frame(state=c(1,2),prior=c(0.5,0.5),lambda=c(.08,.08),mu=c(0.001,0.001),
                   psi=c(.01,.09))

saasi_result <- saasi(bird.orders,pars,ace_qij)

## ----transition adj-----------------------------------------------------------

adj_qij <- q_adjust(ace_qij,"1",1/9)
saasi_result_with_adjusted_qij <- saasi(bird.orders,pars,adj_qij)


## ----echo = TRUE, results = 'hide', fig.show = 'hide'-------------------------
plot(bird.orders, label.offset = 1)

color <- c("#E41A1C", "#377EB8")

# Map tip states to colors
tip_colors <- color[bird.orders$tip.state]
tiplabels(bg = tip_colors, cex = 1, adj = 1, pch = 21)

nodelabels(pie = saasi_result, 
           piecol = color[1:ncol(saasi_result)], 
           cex = 0.2)# plotting


## ----echo = TRUE, results = 'hide', fig.show = 'hide'-------------------------
plot(bird.orders, label.offset = 1)

color <- c("#E41A1C", "#377EB8")

# Map tip states to colors
tip_colors <- color[bird.orders$tip.state]
tiplabels(bg = tip_colors, cex = 1, adj = 1, pch = 21)

nodelabels(pie = saasi_result_with_adjusted_qij, 
           piecol = color[1:ncol(saasi_result_with_adjusted_qij)], 
           cex = 0.2)# plotting


