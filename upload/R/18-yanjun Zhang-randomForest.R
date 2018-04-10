#将数据集分为训练集和测试集,并查看数据集基本属性。数据为R自带IRIS数据
ind<-sample(2,nrow(iris),replace=TRUE,prob=c(0.7,0.3))
set.seed(100)
train<-iris[ind==1,]
test<-iris[ind==2,]
str(train)
str(test)

#选取randomforest –mtry节点值，对应误差最小为2，一般可默认。通常也是2记得
#mtry指定节点中用于二叉树的变量个数，默认情况下数据集变量个数的二次方根（分类模型）或三分之一（预测模型)
library(randomForest)
n <- length(names(train))
set.seed(100)
for (i in 1:(n-1)){
  mtry_fit <- randomForest(Species~.,data=train,mtry=i)
  err <- mean(mtry_fit$err.rate)
  print(err)
}

#之后选择ntree值，ntree指定随机森林所包含的决策树数目，默认为500；.在600左右时，模型内误差基本稳定，故取ntree=600
set.seed(100)
ntree_fit <- randomForest(Species~.,data=train,mtry=2,ntree=1000)
plot(ntree_fit)

#看结果
set.seed(100)
rf <- randomForest(Species~.,data=train,mtry=2,ntree=600,importance=TRUE)
rf

#看重要性
importance <-importance(x=rf)
importance
set.seed(100)
varImpPlot(rf)

#最后验证并预测
pred1 <-predict(rf,data=train)
freq1 <- table(pred1,train$Species)

#验证矩阵中迹占整体情况
sum(diag(freq1))/sum(freq1)

plot(margin(rf,test$Species))
