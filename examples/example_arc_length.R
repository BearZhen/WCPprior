library(sp)
theta <- seq(from=pi,to=3*pi/2,length.out = 50)
line <- Line(cbind(sin(theta),1+ cos(theta)))
sline <- sp::SpatialLines(list(Lines(list(line),ID="1")))

plot(sline)

coords <- sline@lines[[1]]@Lines[[1]]@coords


library(PCWprior)
# Needs a matrix nx2
compute_partial_arc_lengths(coords)

arclengths <- compute_partial_arc_lengths(coords)[,3]

true_arclengths <- (theta-pi)

print("Difference between the true and the piecewiselinear approximation")
sum((true_arclengths - arclengths)^2)
