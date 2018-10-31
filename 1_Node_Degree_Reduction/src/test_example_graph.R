g1 <- make_ring(5)
V(g1)$name <- letters[1:5]

g2 <- make_ring(5)
V(g2)$name <- letters[6:10]

g3 <- make_ring(5)
V(g3)$name <- letters[11:15]

g4 <- make_tree(10, mode = "undirected")
V(g4)$name <- LETTERS[1:10]

g5 <- make_tree(15, mode = "undirected")
V(g5)$name <- LETTERS[11:25]

g <- g1 + g2 + g3 + g4# + g4 # + path("a", "m", "o", "f", "l")  # + edge("a", "A") + edge("c", "C")
plot(g)

plot(GetLCC(g))

graph <- g
components <- components(g)
components
groups(clu)
which.max(components$csize)

plot(lcc)
