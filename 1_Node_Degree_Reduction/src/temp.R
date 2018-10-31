
if (!require("pacman")) install.packages("pacman")
pacman::p_load_gh("trinker/qdapTools")

g1 <- new.env(hash = TRUE)
g2 <- new.env(hash = TRUE)

init1 <- function(s) { g1[[s]] <<- 0L }
init2 <- function(s) { g2[[s]] <<- 0L }

count1 <- function(s) { g1[[s]] <<- g1[[s]] <<- g1[[s]] + 1L }
count2 <- function(s) { g2[[s]] <<- g2[[s]] <<- g2[[s]] + 1L }

lapply(as_ids(V(mm.graph)), init1)
lapply(as_ids(V(mm.graph)), init2)

first <- head(as_ids(V(mm.graph)))
lapply( X = first, FUN = count1)
lapply( X = as_ids(V(mm.graph)), FUN = count1)
lapply( X = first, FUN = count1)

lapply( X = first, FUN = count2)

for(v in as_ids(V(mm.graph))) {
  if(g1[[v]] > 1)
    cat(v, ": ", g1[[v]], "\n")
}

list1 <- unlist(as.list(g1))
df1 <- data.frame(key = names(list1), value = list1, row.names = NULL)
df1$Group <- "1"

list2 <- unlist(as.list(g2))
df2 <- data.frame(key = names(list2), value = list2, row.names = NULL)
df2$Group <- "2"

df <- rbind(df1, df2)

GetSubcomponents(graph = mm.graph, )
