set.seed(100)
data("sonar_mines")
x = sonar_mines[,-1]
y = sonar_mines$y
class1 = "R"
class2 ="M"

y = as.character(y)
y[y==class1]=1
y[y==class2]=0
y = as.numeric(y)
y = matrix(y,ncol = 1)

data_test_x = x
data_test_x = as.matrix(data_test_x)
data_test_y = as.matrix(y)
data_test_x = scale(data_test_x)
PLStests(data_test_y,data_test_x,family="binomial")

## add power2
x_x_2 = cbind(data_test_x,data_test_x^2)
x_x_2 = scale(x_x_2)
PLStests(data_test_y,x_x_2,family="binomial")
