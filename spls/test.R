

test1 = c(0,1,2,3)
test2 = cbind(c(2,3,4,5), c(3,4,5,6))
print(test2)
test2 = test2[rowSums(test2 > 2) > 1,]
print(test2)
