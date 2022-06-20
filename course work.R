library(readxl)
library(copula) 
library(openxlsx)
data=read_excel("~/Desktop/data_10shares.xlsx")
data=data*100 #переход к %
num_of_shares=ncol(data) #количество цб
names_shares=c("Аэрофлот","Северсталь","Газпром","Лукойл","Магнит", "МТС", "НЛМК", "Роснефть", "Сбербанка", "ВТБ")

#Генерация портфелей
G=10000 
X=matrix(ncol=num_of_shares,nrow=G) 
for (k in 1:G) {
  p=c(0,1,runif(n=num_of_shares-1,min=0,max=1)) 
  q=t(matrix(sort(p))) 
  n=ncol(q) 
  w=c(rep(0,times=n-1)) 
  for (j in 1:n-1) { 
    w[j]=q[j+1]-q[j]
  }
  X[k,]=w 
}

#Моделирование методом Марковица
er_markovits=matrix(nrow=1, ncol=num_of_shares) 
for (i in 1:num_of_shares){
  er_markovits[i]=mean(data[ ,i])
}

std_markovits=matrix(nrow=1,ncol=G)
portfolio_r_markovits=matrix(nrow=1,ncol=G)
covarinace_markovits=cov(data)

for (j in 1:G){
  portfolio_r_markovits[j]=sum(er_markovits*X[j,])
  risk_nakop=0
  for (k in 1:num_of_shares){
    for (i in 1:num_of_shares){
      risk_nakop=risk_nakop+X[j,k]*X[j,i]*covarinace_markovits[k,i]
    }
  }
  std_markovits[j]=sqrt(risk_nakop)
}

x11()
plot(std_markovits, portfolio_r_markovits*100, xlab="Портфельный риск",ylab="Ожидаемые доходности", main="Достижимое множество Марковица") 

min_std_markovits=order(std_markovits, decreasing = FALSE)
portfolio_min_std_markovits=X[min_std_markovits[1], ]

rsigma_markovits=portfolio_r_markovits/std_markovits
max_rsigma_markovits=order(rsigma_markovits, decreasing = TRUE)
portfolio_max_rsigma_markovits=X[max_rsigma_markovits[1], ] 

print(portfolio_min_std_markovits)
print(portfolio_max_rsigma_markovits)

#Создание копул

cdf_func=pobs(data)

norm_copula=normalCopula(dim=num_of_shares,dispstr="ex") #Копула Гаусса
t_copula=tCopula(dim=num_of_shares,dispstr="ex") #Kопула Стьюдента
gumbel_copula=gumbelCopula(dim=num_of_shares,param=2) #Kопула Гумбеля
clayton_copula=claytonCopula(dim=num_of_shares,param=2) #Kопула Клейтона
#Подгонка
norm_fit=fitCopula(cdf_func,copula=norm_copula)
t_fit=fitCopula(cdf_func,copula=t_copula)
gumbel_fit=fitCopula(cdf_func,copula=gumbel_copula)
clayton_fit=fitCopula(cdf_func,copula=clayton_copula)
#Выбор оптимальной копулы
norm_fit@loglik
t_fit@loglik
gumbel_fit@loglik
clayton_fit@loglik

#Моделирование доходностей с многомерным распределением на основе лучшей копулы
Z=1000 #количество симуляций
r.week=matrix(nrow=Z,ncol=G) 
r.month=matrix(nrow=Z,ncol=G) 
r_c.week=matrix(nrow=Z,ncol=num_of_shares) 
r_c.month=matrix(nrow=Z,ncol=num_of_shares)
for (l in 1:Z) {
  D=30 
  copula_gen=rCopula(n=D,copula=t_fit@copula) 
  r_matrix=matrix(nrow=D,ncol=num_of_shares) 
  r_matrix_c.week=matrix(nrow=1,ncol=num_of_shares)
  r_matrix_c.month=matrix(nrow=1,ncol=num_of_shares)
  for (j in 1:num_of_shares) {
    fact_cdf=ecdf(data[,j]) 
    r_matrix[,j]=quantile(fact_cdf,t(copula_gen[,j])) 
    r_matrix_c.week[j]=(prod(1+r_matrix[1:7,j]/100)-1)*100 
    r_matrix_c.month[j]=(prod(1+r_matrix[1:30,j]/100)-1)*100 
  }
  #Считаем CVaR и ожидаемую доходность
  r_portfolio=matrix(ncol=G,nrow=D) 
  for (i in 1:G) {
    ZZ=0
    for (j in 1:num_of_shares){
      ZZ=ZZ+r_matrix[,j]*X[i,j] 
      r_portfolio[,i]=ZZ
    }
  }
  r_portfolio_week=matrix(ncol=G,nrow=1)
  r_portfolio_month=matrix(ncol=G,nrow=1)
  for (i in 1:G) {
    r_portfolio_week[,i]=(prod((r_portfolio[1:7,i]/100)+1)-1)*100
    r_portfolio_month[,i]=(prod((r_portfolio[1:30,i]/100)+1)-1)*100
  }
  r.week[l,]=r_portfolio_week
  r.month[l,]=r_portfolio_month
  r_c.week[l,]=r_matrix_c.week
  r_c.month[l,]=r_matrix_c.month
}

mean_r_week=matrix(ncol=G,nrow=1)
mean_r_month=matrix(ncol=G,nrow=1)
CVaR_week=matrix(ncol=G,nrow=1)
CVaR_month=matrix(ncol=G,nrow=1)
VaR_week=matrix(ncol=G,nrow=1)
VaR_month=matrix(ncol=G,nrow=1)
gamma_=0.95 
VaR_bound=(1-gamma_)*D
for (i in 1:G) {
  mean_r_week[,i]=mean(r.week[,i]) 
  mean_r_month[,i]=mean(r.month[,i]) 
  sort_r_week=sort(r.week[,i]*(-1),decreasing=TRUE)
  sort_r_month=sort(r.month[,i]*(-1),decreasing=TRUE)
  VaR_week[i]=sort_r_week[VaR_bound]
  VaR_month[i]=sort_r_month[VaR_bound]
  CVaR_week[i]=mean(sort_r_week[1:VaR_bound])
  CVaR_month[i]=mean(sort_r_month[1:VaR_bound])
}
a.par=par(mfrow=c(1,2))
plot(CVaR_week, mean_r_week, xlab="CVaR",ylab="Ожидаемая доходность", main="Достижимое множество недельные данные")
plot(CVaR_month, mean_r_month, xlab="CVaR",ylab="Ожидаемая доходность", main="Достижимое множество месячные данные")

#Поиск оптимального портфеля по максимуму r/CVaR
coeff_week=mean_r_week/CVaR_week
max_coeff_week=order(coeff_week,decreasing = TRUE)
max_coeff_week_portfolio=X[max_coeff_week[1],] 
char_max_coeff_week=c(mean_r_week[max_coeff_week[1]], CVaR_week [max_coeff_week[1]], VaR_week[max_coeff_week[1]], coeff_week[max_coeff_week[1]])
names(char_max_coeff_week)=c("Доходность", "CVaR ","VaR ","Коэффициент")

coeff_month=mean_r_month/CVaR_month
max_coeff_month=order(coeff_month,decreasing = TRUE)
max_coeff_month_portfolio=X[max_coeff_month[1],] 
char_max_coeff_month=c(mean_r_month[max_coeff_month[1]], CVaR_month[max_coeff_month[1]], VaR_month[max_coeff_month[1]], coeff_month[max_coeff_month[1]])
names(char_max_coeff_month)=c("Доходность", "CVaR ","VaR ","Коэффициент")

# Портфель минимального риска по CVaR
min_CVaR_week=order(CVaR_week, decreasing = FALSE)
min_CVaR_week_portfolio=X[min_CVaR_week[1],]
char_min_CVaR_week=c(mean_r_week[min_CVaR_week[1]],CVaR_week[min_CVaR_week[1]],VaR_week[min_CVaR_week[1]], coeff_week[min_CVaR_week[1]])
names(char_min_CVaR_week)=c("Доходность", "CVaR ","VaR ","Коэффициент")

min_CVaR_month=order(CVaR_month, decreasing = FALSE)
min_CVaR_month_portfolio=X[min_CVaR_month[1],]
char_min_CVaR_month=c(mean_r_month[min_CVaR_month[1]],CVaR_month[min_CVaR_month[1]],VaR_month[min_CVaR_month[1]], coeff_month[min_CVaR_month[1]])
names(char_min_CVaR_month)=c("Доходность", "CVaR ","VaR ","Коэффициент")

print(min_CVaR_week_portfolio)
print(max_coeff_week_portfolio)
print(min_CVaR_month_portfolio)
print(max_coeff_month_portfolio)

std_markovits_2=std_markovits
for (i in 1:G){
  std_markovits_2[i]=std_markovits_2[i]+0.5*(sum(X[i,]^2))
}
x11()
plot(std_markovits_2, portfolio_r_markovits*100, xlab="Портфельный риск",ylab="Ожидаемые доходности", main="Достижимое множество Марковица") 

min_std_markovits_2=order(std_markovits_2, decreasing = FALSE)
portfolio_min_std_markovits_2=X[min_std_markovits_2[1], ]
portfolio_min_std_markovits_2

CVaR_week_2=CVaR_week
CVaR_month_2=CVaR_month

for (i in 1:G){
  CVaR_week_2[i]=CVaR_week_2[i]+0.5*(sum(X[i,]^2))
}

for (i in 1:G){
  CVaR_month_2[i]=CVaR_month_2[i]+0.5*(sum(X[i,]^2))
}

a.par=par(mfrow=c(1,2))
plot(CVaR_week_2, mean_r_week, xlab="CVaR",ylab="Ожидаемая доходность", main="Достижимое множество недельные данные")
plot(CVaR_month_2, mean_r_month, xlab="CVaR",ylab="Ожидаемая доходность", main="Достижимое множество месячные данные")

min_CVaR_week_2=order(CVaR_week_2, decreasing = FALSE)
min_CVaR_week_portfolio_2=X[min_CVaR_week_2[1],]
min_CVaR_week_portfolio_2

min_CVaR_month_2=order(CVaR_month_2, decreasing = FALSE)
min_CVaR_month_portfolio_2=X[min_CVaR_month_2[1],]
min_CVaR_month_portfolio_2

