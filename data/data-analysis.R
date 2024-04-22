source("lin-anc-ts.R")

require(MASS) # for geyser data
# analyze geyser data
la_geyser <- lin.anc.ts(geyser, 6)
instant.p.val(la_geyser)
instant.graph(la_geyser)
summary.p.val(la_geyser)
summary.graph(la_geyser)

# shift waiting such that it is waiting after erruption
geyser2 <- data.frame(waiting = geyser$waiting[-1], duration = geyser$duration[-nrow(geyser)])
la_geyser2 <- lin.anc.ts(geyser2, 6)
instant.p.val(la_geyser2)
instant.graph(la_geyser2)
summary.p.val(la_geyser2)
summary.graph(la_geyser2)


# analyze gas furnace data
# data should be stored in data/gas-furnace.csv available from https://openmv.net/info/gas-furnace
gas_furnace <- read.csv("data/gas-furnace.csv")

la_gas <- lin.anc.ts(gas_furnace, 6)
instant.p.val(la_gas)
instant.graph(la_gas)
summary.p.val(la_gas)
summary.graph(la_gas)

# analyze dairy data
# data should be stored in data/data_butter.txt and data/data_cheddar.txt
# create data frame
milk <- data.frame(butter = read.table("data/data_butter.txt")$V1, cheddar = read.table("data/data_cheddar.txt")$V1)

la_milk <- lin.anc.ts(milk, 6)
instant.p.val(la_milk)
instant.graph(la_milk)
summary.p.val(la_milk)
summary.graph(la_milk)