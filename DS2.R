install.packages("UsingR")
library(UsingR)
data("father.son",package="UsingR")
head(father.son)
mean(father.son$sheight)
tapply(father.son$sheight, round(father.son$fheight), mean)

X = matrix(1:1000,100,10)
X[25,3]
x=1:10
A = cbind(x,x2=2*x,x3=3*x,x4=4*x,x5=5*x)
sum(A[7,])

X = matrix(c(3,2,1,5,4,2,-1,0,-5,2,5,0,1,-1,-5,1),4,4)
y <- matrix(c(10,5,7,4),4,1)
solve(X)%*%y

a <- matrix(1:12, nrow=4)
b <- matrix(1:15, nrow=3)
a
b
a%*%b
sum(a[3,]*b[,2])

X <- matrix(c(1,1,1,1,0,0,1,1),nrow=4)
rownames(X) <- c("a","a","b","b")
beta <- c(5, 2)
fitted = X %*% beta
fitted[ 1:2, ]
X <- matrix(c(1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,1,1),nrow=6)
rownames(X) <- c("a","a","b","b","c","c")
beta <- c(10,3,-3)
fitted = X %*% beta
 
g = 9.8
h0 = 56.67
v0 = 0
n = 25
tt = seq(0,3.4,len=n) 
set.seed(1)
ME = replicate(100000,{y = h0 + v0 *tt - 0.5* g*tt^2 + rnorm(n,sd=1)
    X = cbind(1,tt,tt^2)
    A = solve(crossprod(X))%*%t(X)
    -2 * (A%*%y)[3]
})
sd(ME)

library(UsingR)
x = father.son$fheight
y = father.son$sheight
n = length(y)
N =  50
B = 10000
set.seed(1)
SE = replicate(B,{index = sample(n,N)
    sampledat = father.son[index,]
    x = sampledat$fheight
    y = sampledat$sheight
    betahat = lm(y~x)$coef[2]
    return(betahat)
    })
sd(SE)
mean((y - mean(y))*(x-mean(x)))

library(UsingR)
x = father.son$fheight
y = father.son$sheight
n = length(y)
N = 50
set.seed(1)
index = sample(n,N)
sampledat = father.son[index,]
x = sampledat$fheight
y = sampledat$sheight
betahat = lm(y~x)$coef[2]
fit = lm(y~x)
yf=fit$fitted.values
SSE = sum((y-yf)^2)
SSE
sigma2 = SSE/48
X = cbind(1,x)
X2=solve(t(X)%*%X)
sqrt(diag(X2)*sigma2)

nx = 5
ny = 7
X = cbind(rep(1,nx + ny),rep(c(0,1),c(nx, ny)))
(t(X)%*%X)[1,1]

install.packages("contrast")
library(contrast)
species <- factor(c("A","A","B","B"))
condition <- factor(c("control","treated","control","treated"))
model.matrix(~ species + condition)
y = rnorm(4)
fit = lm(y ~ species + condition)
contrast(fit, list(species="B",condition="control"), list(species="A",condition="treated"))$X

url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/spider_wolff_gorb_2013.csv"
filename <- "spider_wolff_gorb_2013.csv"
library(downloader)
if (!file.exists(filename)) download(url, filename)
spider <- read.csv(filename, skip=1)
fit = lm(friction ~ leg+type, data=spider)
contrast(fit,list(leg="L4",type="pull"), list(leg="L2",type="pull"))

X <- model.matrix(~ type + leg, data=spider)
(Sigma <- sum(fit$residuals^2)/(nrow(X) - ncol(X)) * solve(t(X) %*% X))
C <- matrix(c(0,0,-1,0,1),1,5)
sqrt((0.0011819981+0.0020871318)-(2*0.0006389179))
sqrt(C %*% Sigma %*% t(C))

spider$log2friction <- log2(spider$friction)
boxplot(log2friction ~ type*leg, data=spider)

fit = lm(log2friction~type*leg, data=spider)
summary(fit)
anova(fit)
contrast(fit, list(type="pull",leg="L2"), list(type="pull",leg="L1"))
contrast(fit, list(type="push",leg="L2"), list(type="push",leg="L1"))

N <- 40
p <- 4
set.seed(1)
fval = replicate(1000,{
  group <- factor(rep(1:p,each=N/p))
  X <- model.matrix(~ group)
  Y <- rnorm(N,mean=42,7)
  mu0 <- mean(Y)
  initial.ss <- sum((Y - mu0)^2)
  s <- split(Y, group)
  after.group.ss <- sum(sapply(s, function(x) sum((x - mean(x))^2)))
  (group.ss <- initial.ss - after.group.ss)
  group.ms <- group.ss / (p - 1)
  after.group.ms <- after.group.ss / (N - p)
  group.ms / after.group.ms
})
mean(fval)
hist(fval, col="grey", border="white", breaks=50, freq=FALSE)
df1=p - 1
df2=N - p
xs <- seq(from=0,to=6,length=100)
lines(xs, df(xs, df1 = p - 1, df2 = N - p), col="red")

sex <- factor(rep(c("female","male"),each=4))
trt <- factor(c("A","A","B","B","C","C","D","D"))
X <- model.matrix( ~ sex + trt)
X
qr(X)$rank
Y <- 1:8
makeYstar <- function(a,b) Y - X[,2] * a - X[,5] * b
fitTheRest <- function(a,b) {
  Ystar <- makeYstar(a,b)
  Xrest <- X[,-c(2,5)]
  betarest <- solve(t(Xrest) %*% Xrest) %*% t(Xrest) %*% Ystar
  residuals <- Ystar - Xrest %*% betarest
  sum(residuals^2)
}
fitTheRest(1,2)
expand.grid(1:3,1:3)
betas = expand.grid(-2:8,-2:8)
rss = apply(betas,1,function(x) fitTheRest(x[1],x[2]))
rss
min(rss)
fitTheRest(1,5)
library(rafalib)
themin=min(rss)
plot(betas[which(rss==themin),])

fit <- lm(friction ~ type + leg, data=spider)
betahat <- coef(fit)
Y <- matrix(spider$friction, ncol=1)
X <- model.matrix(~ type + leg, data=spider)
QR = qr(X)
Q = qr.Q(QR)
head(Q)
R = qr.R(QR)
head(R)
crossprod(Q,Y)
