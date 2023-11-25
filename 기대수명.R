# 파일 경로 설정 #
# 아래에도 파일 경로를 적절하게 변경해야 정상적으로 실행이 가능합니다. #
filePath = ""

# 단순선형회귀 함수 #
simpleLinear = function(fileName, xName) {
  library(bbmle)
  temp = paste0(filePath, fileName, ".csv")
  data = read.csv(file = temp, header = T)
  # 상관계수 #
  cor(data[,2],data$기대수명)
  ## 단순회귀 ##
  # 최소제곱추정량 #
  result = lsfit(x=data[,2],y=data$기대수명)
  ls.print(result)
  
  # 그래프 그리기 #
  plot(x=data[,2],y=data$기대수명,xlab = as.character(xName), ylab = "Life Expectancy")
  abline(result)
}

simpleLinear("영아사망률","Infant Death Rate") # 0.8268
simpleLinear("인터넷","Internet") # 0.7162
simpleLinear("GNI","GNI per Capita") # 0.4255
abline(v=1e+4)
simpleLinear("GNI_low","GNI_low per Capita") # 0.4143
simpleLinear("GNI_high","GNI_high per Capita") # 0.2827

simpleLinear("AIDS","AIDS") # 0.2652 # 이분산
simpleLinear("AIDS2","AIDS2") # 0.2574 # 자유도 추가
simpleLinear("CO2","CO2 per Capita") # 0.2909 # GNI랑 비슷함
simpleLinear("흡연율","Cigarette") # 0.111
simpleLinear("알콜","Alcohol Consumption per Capita") # 0.0637

simpleLinear("재난피해자","Disaster Victims") # 0.069 # 이분산 발생 

simpleLinear("보건지출비","Health Service expenditure") # 0.0158 # 아래와 반대의 결과
simpleLinear("보건지출비2","Health Service expenditure2") # 0.183 # 극단적인 값 제외

simpleLinear("산림면적비율","Forest Area ratio") # 0

### 1인당GNI ### (최우추정량)
temp = paste0(filePath, "GNI.csv")
data = read.csv(file = temp, header = T)
#install.packages("bbmle")
library(bbmle)
y = data$기대수명
x = data$GNI
n = length(y)
fn = function(b0, b1, sigma) { 
  (n/2) * log(sigma^2) + 1 / (2*sigma^2) * (sum ((y - b0 - b1 * x) ^ 2)) 
}
res = mle2(fn, start = list(b0 = 68, b1 = 0, sigma = 0.1), data = data)
summary(res)

### CO2 ### (최우추정량)
temp = paste0(filePath, "CO2.csv")
data = read.csv(file = temp, header = T)
#install.packages("bbmle")
library(bbmle)
y = data$기대수명
x = data$CO2
n = length(y)
fn = function(b0, b1, sigma) { 
  (n/2) * log(sigma^2) + 1 / (2*sigma^2) * (sum ((y - b0 - b1 * x) ^ 2)) 
}
res = mle2(fn, start = list(b0 = 68, b1 = 0, sigma = 1), data = data)
summary(res)

predicted <- coef(res)["sigma"] + coef(res)["b0"] * x
# 잔차 계산
residuals <- y - predicted
# 잔차 히스토그램 그리기
hist(residuals)

### 다중 공산성 탐지 ### (GNI, CO2)
temp = paste0(filePath, "GNI+CO2.csv")
GC = read.csv(file = temp, header = T,row.names = 1)
GC.res = lm(기대수명~GNI+CO2, data = GC)
summary(GC.res)
cor(GC)
library(HH)
vif(GC[,c(1,2)])

### 이분산 GLS ### (AIDS)S
temp = paste0(filePath, "AIDS2.csv")
AIDS2 = read.csv(file = temp, header = T)
AIDS2.res = lm(기대수명~AIDS,data = AIDS2)
residuals = resid(AIDS2.res)
yhat = predict(AIDS2.res, interval = "none")
plot(x=yhat, y=residuals)

library(nlme)
AIDS2.gls.res = gls(model = 기대수명~AIDS, data = AIDS2, weights = varFixed(~기대수명))
AIDS2.gls.res


### 이분산 GLS ### (CO2)
temp = paste0(filePath, "CO2.csv")
CO2.df = read.csv(file = temp, header = T)
CO2.res = lm(기대수명~CO2,data = CO2.df)
residuals = resid(CO2.res)
yhat = predict(CO2.res, interval = "none")
plot(x=yhat, y=residuals)

library(nlme)
CO2.gls.res = gls(model = 기대수명~CO2, data = CO2.df, weights = varFixed(~기대수명))
summary(CO2.gls.res)
CO2.gls.res
#simpleLinear("CO2","AIDS2") # 0.2574 # 자유도 추가

### 이분산 GLS ### (재난피해자)
temp = paste0(filePath, "재난피해자.csv")
dis.df = read.csv(file = temp, header = T)
dis.res = lm(기대수명~재난피해자,data = dis.df)
residuals = resid(dis.res)
yhat = predict(dis.res, interval = "none")
plot(x=yhat, y=residuals)

library(nlme)
dis.gls.res = gls(model = 기대수명~재난피해자, data = dis.df, weights = varFixed(~기대수명))
dis.gls.res

### 다중회귀 ### (GNI, 인터넷, CO2)
temp = paste0(filePath, "multiple3.csv")
multiple3 = read.csv(file = temp, header = T)
ls3.res = lm(formula = 기대수명~GNI+인터넷+CO2,data = multiple3)
summary(ls3.res)
BIC(ls3.res)
### 다중회귀 ### (GNI, 인터넷)
temp = paste0(filePath, "multiple3.csv")
multiple2 = read.csv(file = temp, header = T)
ls2.res = lm(formula = 기대수명~GNI+인터넷,data = multiple3)
summary(ls2.res)
BIC(ls2.res)

### 다중회귀 최우추정 ### (GNI, 인터넷, CO2)
library(bbmle)
y = multiple3$기대수명
x1 = multiple3$GNI
x2 = multiple3$인터넷
x3 = multiple3$CO2
n = length(y)
fn3 = function(b0, b1, b2, b3, sigma) { 
  (n/2) * log(sigma^2) + 1 / (2*sigma^2) * (sum ((y - b0 - b1 * x1 - b2 * x2 - b3 * x3) ^ 2))
}
res3 = mle2(fn3 ,start = list(b0 = 60, b1 = 0, b2 = 0, b3 = 0, sigma = 1))
summary(res3)

### 다중회귀 ### (GNI, 인터넷, CO2, 흡연율)
temp = paste0(filePath, "multiple4.csv")
multiple4 = read.csv(file = temp, header = T)
ls4.res = lm(formula = 기대수명~GNI+인터넷+CO2+흡연율,data = multiple4)
summary(ls4.res)
BIC(ls4.res)

### 다중회귀 ### (GNI, 인터넷, 흡연율)
temp = paste0(filePath, "multiple4.csv")
multiple4 = read.csv(file = temp, header = T)
ls41.res = lm(formula = 기대수명~GNI+인터넷+흡연율,data = multiple4)
summary(ls41.res)
BIC(ls41.res)

### 다중회귀 최우추정 ### (GNI, 인터넷, CO2,흡연율)
library(bbmle)
temp = paste0(filePath, "multiple4.csv")
multiple4 = read.csv(file = temp, header = T)
y = multiple4$기대수명
x1 = multiple4$GNI
x2 = multiple4$인터넷
x3 = multiple4$CO2
x4 = multiple4$흡연율
n = length(y)
fn4 = function(b0, b1, b2, b3, b4, sigma) { 
  (n/2) * log(sigma^2) + 1 / (2*sigma^2) * (sum ((y - b0 - b1 * x1 - b2 * x2 - b3 * x3 - b4*x4) ^ 2))
}
res4 = mle2(fn4 ,start = list(b0 = 60, b1 = 0, b2 = 0, b3 = 0, b4 = 0, sigma = 1))
summary(res4)

### 다중회귀 최우추정 ### (GNI, 인터넷, 흡연율)
library(bbmle)
temp = paste0(filePath, "multiple4.csv")
multiple4 = read.csv(file = temp, header = T)
y = multiple4$기대수명
x1 = multiple4$GNI
x2 = multiple4$인터넷
x3 = multiple4$흡연율
n = length(y)
fn41 = function(b0, b1, b2, b3, sigma) { 
  (n/2) * log(sigma^2) + 1 / (2*sigma^2) * (sum ((y - b0 - b1 * x1 - b2 * x2 - b3*x3) ^ 2))
}
res41 = mle2(fn41 ,start = list(b0 = 60, b1 = 0, b2 = 0, b3 = 0, sigma = 1))
summary(res41)


# """폐기"""
### 다중회귀 ### (GNI, 인터넷, 흡연율, AIDS)
temp = paste0(filePath, "multiple_A.csv")
multiple_A = read.csv(file = temp, header = T)
lsA.res = lm(formula = 기대수명~GNI+인터넷+흡연율+AIDS,data = multiple_A)
summary(lsA.res)
BIC(lsA.res)

### 다중회귀 ### (GNI, 인터넷, 흡연율, 보건지출비)
temp = paste0(filePath, "multiple5.csv")
multiple5 = read.csv(file = temp, header = T)
ls5.res = lm(formula = 기대수명~GNI+인터넷+흡연율+보건지출비,data = multiple5)
summary(ls5.res)
BIC(ls5.res)
