require(MASS)
source("~/Documents/ETH/PhD/SVAR-Ancestor/lin-anc-ts.R")

# shift waiting such that it is waiting after erruption
geyser2 <- data.frame(waiting = geyser$waiting[-1], duration = geyser$duration[-nrow(geyser)])

la_geyser <- lin.anc.ts(geyser, 6)
la_geyser2 <- lin.anc.ts(geyser2, 6)

summary.graph(la_geyser)
summary.graph(la_geyser2)
# same summary graph

instant.graph(la_geyser)
instant.graph(la_geyser2)
# only latter shows instant effect


milk <- data.frame(butter = read.table("data/data_butter.txt")$V1, cheddar = read.table("data/data_cheddar.txt")$V1)

la_milk <- lin.anc.ts(milk, 6)
instant.graph(la_milk)
summary.graph(la_milk)

gas_furnace <- read.csv("data/gas-furnace.csv")

la_gas <- lin.anc.ts(gas_furnace, 6)
instant.graph(la_gas)
summary.graph(la_gas)
