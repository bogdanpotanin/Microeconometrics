# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Урок 10. Модель Тобина и усеченная регрессия
# --------

# Отключим scientific notation
options(scipen = 999)

#---------------------------------------------------
# Часть 1. Оценивание параметров
#---------------------------------------------------

# Подключим дополнительные библиотеки
library("mvtnorm")                       # симуляции из многомерного
                                         # нормального распределения

library("numDeriv")                      # численное дифференцирование

library("stringr")                       # работа со строками
library("tm")

library("AER")                           # tobit модель и тест на overdispersion 
library("VGAM")                          # Модель tobit, второй метод
library("crch")                          # Модель tobit, третий метод
library("truncreg")                      # Регрессия с усечением

library("hpa")                           # моменты усеченного 
                                         # нормального распределения

#---------------------------------------------------
# Часть 1. Оценивание параметров
#---------------------------------------------------

#*******************************
# Пользовательский уровень
#*******************************

set.seed(123)
# Зададим количество наблюдений
n <- 10000
# Истинные значения регрессионных коэффициентов
beta_0 <- 6
beta_1 <- 2
beta_2 <- 1.5
beta_3 <- 1
beta <- c(beta_0, beta_1, beta_2, beta_3)
# Создадим независимые переменные из 
# многомерного нормального распределения
X <- rmvnorm(n,                                 # симулируем n наблюдений из многомерного
                                                # нормального распределения
             c(0, 0, 0),                        # с нулевым вектором математических ожиданий и
             matrix(c(1, 0.2, 0.3,              # следующей ковариационной матрице
                      0.2, 1, -0.1,
                      0.3, -0.1, 1),
                    ncol = 3,
                    byrow = FALSE))
# Случайная ошибка
sigma = 2.5
epsilon <- rnorm(n, mean = 0, sd = sigma)
# Не усеченная зависимая переменная
y_star <- beta_0 + beta_1 * X[, 1] + beta_2 * X[, 2] + 
          beta_3 * X[, 3] + epsilon
# Точки усечения
tr_left <- 2
tr_right <- 12
# Усеченная зависимая переменная
y <- y_star
y[y_star <= tr_left] <- tr_left
y[y_star >= tr_right] <- tr_right
# Сюжет:
# Представим, что вы изучаете факторы, влияющие на объем
# выдаваемого индивиду кредита. Представим, что вследствие законодательно 
# введенных ограничений банк выдает лишь кредиты объемом от 2 до 12 условных 
# единиц. В таком случаи переменная y_star отражает то, сколько банк был бы 
# готов дать индивиду, если бы ограничения отсутствовали, а y - то, сколько
# дает фактически.
h <- data.frame("credit" = y,                        # объем выданного кредита
                "age" = X[, 1],                      # возраст
                "income" = X[, 2],                   # доход
                "story" = X[, 3])                    # кредитная история
# Оценим параметры модели
  # при помощи модели Тобина
model_tobit <- crch(credit ~ age + income + story,   # формула
                    data = h,                        # данные
                    left = tr_left,                  # нижнее (левое) усечение
                    right = tr_right)                # верхнее (правое) усечение
summary(model_tobit)                                 # посмотрим результат
est_tobit <- coef(model_tobit)                       # достанем оценки
coef_tobit <- est_tobit[-length(est_tobit)]          # оценки регрессионных
                                                     # коэффициентов
sigma_tobit <- exp(est_tobit[length(est_tobit)])     # достаем оценку 
                                                     # стандартного отклонения
  # при помощи усеченной регрессии
model_trunc <- crch(credit ~ age + income + story,   # формула
                    data = h[(h$credit > 2) &        # отбираем лишь кредиты
                             (h$credit < 12), ],     # между 2 и 12,                        
                    left = tr_left,                  # нижнее (левое) усечение
                    right = tr_right,                # верхнее (правое) усечение
                    truncated = TRUE)                # усеченная регрессия               
summary(model_trunc)                                 # посмотрим результат
est_trunc <- coef(model_trunc)                       # достанем оценки
coef_trunc <- est_trunc[-length(est_trunc)]          # оценки регрессионных
                                                     # коэффициентов
sigma_trunc <- exp(est_trunc[length(est_trunc)])     # достаем оценку 
                                                     # стандартного отклонения
  # при помощи МНК по всем наблюдениям
model_lm_all <- lm(credit ~ age + income + story,    # формула
                   data = h)                         # данные
summary(model_lm_all)
coef_lm_all <- coef(model_lm_all)
sigma_lm_all <- sigma(model_lm_all)
  # при помощи МНК по не усеченным наблюдениям
model_lm_tr <- lm(credit ~ age + income + story,     # формула
                  data = h[(h$credit > 2) &          # отбираем лишь кредиты
                           (h$credit < 12), ])       # между 2 и 12
summary(model_lm_tr)
coef_lm_tr <- coef(model_lm_tr)
sigma_lm_tr <- sigma(model_lm_tr)
# Сравним полученные результаты
  # оценки регрессионных коэффициентов
data.frame("True" = beta,
           "Tobit" = coef_tobit,
           "Truncated" = coef_trunc,
           "LS.All" = coef_lm_all,
           "LS.Truncated" = coef_lm_tr)
  # оценки стандартного отклонения
  # случайной ошибки
data.frame("True" = sigma,
           "Tobit" = sigma_tobit,
           "Truncated" = sigma_trunc,
           "LS.All" = sigma_lm_all,
           "LS.Truncated" = sigma_lm_tr)

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 1.1.    Посмотрите, как изменится результат, если:
#         1)     убрать цензурирование слева
#         2)     убрать цензурирование справа
#         3)     убрать цензурирование
#         4)     сдвинуть левую и правую границы
#                ближе друг к другу
# 1.2*.   Придумайте и реализуйте собственный
#         пример на симулированных данных

#*******************************
# Академический уровень
#*******************************

# Запишем функцию правдоподобия
lnL <- function(x, y, X, tr_low, tr_up)
{
  x_n <- length(x)                                       # число оцениваемых
                                                         # параметров
  sigma <- x[x_n]              
  beta <- matrix(x[-x_n], ncol = 1)
  
  X <- cbind(1, X)                                       # добавляем константу
  
  XB <- X %*% beta                                       # считаем произвдение
                                                         # регрессоров на 
                                                         # коэффициенты
  
  n <- nrow(X)
  
  L_vector <- rep(NA, n)
  
  cond_1 <- (y == tr_low)                                # различные условия
  cond_2 <- (y == tr_up)                                 # дают различные вклады
  cond_3 <- ((y > tr_low) & (y < tr_up))                 # в функцию правдоподобия
  
  L_vector[cond_1] <- pnorm(y[cond_1] - XB[cond_1,],     # считаем вклады в
                            0, sigma)                    # функцию правдоподобия
  L_vector[cond_2] <- 1 - pnorm(y[cond_2] - XB[cond_2,], 
                                0, sigma)
  L_vector[cond_3] <- dnorm(y[cond_3] - XB[cond_3,], 0, 
                            sigma)
  
  return(sum(log(L_vector)))                             # возвращаем логарифм
                                                         # функции правдоподобия
}
# Осуществим оптимизацию
x0 <- c(coef_lm_all, sigma_lm_all)                       # в качестве начальной
                                                         # точки возьмем МНК
opt_mle <- optim(par = x0, fn = lnL,                     # запускаем оптимизатор
                 X = X, y = y,
                 tr_low = tr_low, tr_up = tr_up,
                 method = "Nelder-Mead", 
                 control = list("maxit" = 10000,          
                                fnscale = -1,
                                "reltol" = 1e-10))
x1 <- opt_mle$par                                        # сохраняем оценки
# Сравним истинные и предсказанные знаечния
data.frame("beta" = beta, 
           "beta_hat" = x1[-length(x1)])

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 1.1*.   Перепишите функцию правдоподобия, используя
#         репараметризацию Олсена, а именно, замены:
#         beta* = beta / sigma
#         sigma* = 1 / sigma
# 1.2**.  Оцените гипотезы о значимости параметров
#         используя оценки, полученные с использованием
#         репараметризации Олсена

#---------------------------------------------------
# Часть 2. Расчет предсказаний и предельных эффектов
#---------------------------------------------------

#*******************************
# Пользовательский уровень
#*******************************

# Создадим индивида
Boris <- data.frame("age" = 1,
                    "income" = 2,
                    "story" = 3)

# Рассчитаем оценку безусловного математического 
# ожидания зависимой переменной, то есть E(y*)
credit_est <- predict(model_tobit, 
                      newdata = Boris)   

# Рассчитаем оценки вероятностей усечения
  # снизу
prob_est_left <- predict(model_tobit,                    # встроенный метод
                         at = tr_left,                   # нижняя граница усечения
                         type = "probability", 
                         newdata = Boris)                              
prob_est_left <- pnorm(tr_left - credit_est,             # вручную
                       sd = sigma_tobit) 
  # сверху
prob_est_right <- 1 - predict(model_tobit,               # встроенный метод
                              at = tr_right,
                              type = "probability",
                              newdata = Boris)                              
prob_est_right <- 1 - pnorm(tr_right - credit_est,       # вручную
                            sd = sigma_tobit) 

#*******************************
# Продвинутый пользовательский уровень
#*******************************

# Рассчитаем оценку усеченного снизу и сверху 
# математического ожидания зависимой переменной,
# то есть E(y)
epsilon_E <- truncatedNormalMoment(k = 1,                # условное математическое        
               x_lower = tr_left - credit_est,           # ожидание случайной ошибки
               x_upper = tr_right - credit_est,          # E(e|lower < y* < upper)
               mean = 0, sd = sigma_tobit)
epsilon_E <- sigma_tobit ^ 2 *                           # расчет вручную
             (dnorm(tr_left - credit_est, 
                    sd = sigma_tobit) - 
              dnorm(tr_right - credit_est, 
                    sd = sigma_tobit)) / 
             (pnorm(tr_right - credit_est, 
                    sd = sigma_tobit) - 
              pnorm(tr_left - credit_est, 
                    sd = sigma_tobit))
credit_est_tr <- (credit_est + epsilon_E) *              # E(y|lower < y* < upper) *
  (1 - prob_est_left - prob_est_right) +                 # P(lower < y* < upper) +
  prob_est_left * tr_left +                              # P(y > upper) * upper +
  prob_est_right * tr_right                              # P(y < lower) * lower

# Предельный эффект дохода на
ME_credit_income <- coef_tobit["income"]                 # E(y*)
ME_credit_tr_income <- coef_tobit["income"] *            # E(y)
                       (1 - prob_est_right -
                       prob_est_left)
ME_prob_right_income <- coef_tobit["income"] *           # P(y > upper)
                       dnorm(tr_right - credit_est,      
                             sd = sigma_tobit)      
ME_prob_left_income <- -coef_tobit["income"] *           # P(y < lower)
                       dnorm(credit_est - tr_left,      
                             sd = sigma_tobit)

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 1.1.    Для индивида с произвольными характеристиками
#         оцените предельный эффект возраста на:
#         1)     E(y*)
#         2)     P(y > upper)
#         3)     P(y < lower)
#         4)     E(y|lower < y < upper)
#         5)     E(y)
#         6*)    E(y|y > 0)
# 1.2.    Повторите предыдущее задание добавив
#         в модель:
#         1*)    возраст в квадрате
#         2*)    взаимодействие между возрастом
#                и доходом

#---------------------------------------------------
# Часть 3. Гетероскедастичность
#---------------------------------------------------

#*******************************
# Продвинутый пользовательский уровень
#*******************************

# Допустим, что дисперсия случайной ошибки
# зависит от возраста и здоровья
set.seed(123)
h$health <- rnorm(n)                                 # переменная на здоровье
gamma <- c(log(sigma), 0.1, 0.2)                     # коэффициенты уравнения
                                                     # дисперсии
h$credit <- (y_star - epsilon) +                     # учтем гетероскедастичность
             exp(gamma[1] +                          # в зависимой переменной
                 h$age * gamma[2] + 
                 h$health * gamma[3]) * 
            rnorm(n)
h$credit[h$credit > tr_right] <- tr_right
h$credit[h$credit < tr_left] <- tr_left

model_htobit <- crch(credit ~ age + income + story | # основное уравнение
                              age + health,          # уравнение дисперсии
                    data = h,                        # данные
                    left = tr_left,                  # нижнее (левое) усечение
                    right = tr_right,                # верхнее (правое) усечение
                    link.scale = "log")              # тип уравнения дисперсии           
summary(model_htobit) 

# Аргумент dist позволяет оценивать параметры
# модели при альтернативном допущении
# о распределении случайных ошибок
model_htobit_student <- crch(credit ~ 
                     age + income + story |
                     age + health ,         
                     data = h,                        
                     left = tr_left,                  
                     right = tr_right,                
                     link.scale = "log",      
                     dist = "student")               # будем использовать
                                                     # распределение Стьюдента
summary(model_htobit_student) 

# Встроенные функции для расчета предельных эффектов
# и вероятностей работает по аналогии

# Некоторые дополнительные полезные функции
logLik(model_tobit)                                  # правдоподобие
AIC(model_tobit)                                     # AIC
BIC(model_tobit)                                     # BIC

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 1.1.    Сравните модели с нормальным распределением
#         случайных ошибок и распределением Стьюдента
#         по критерию AIC
# 1.2*.   Проверьте, что произойдет, если в уравнение
#         дисперсии будет добавлена незначимая 
#         переменная
# 1.3**.  Для индивида с произвольными характеристиками
#         оцените предельный эффект возраста на:
#         1)     E(y*)
#         2)     P(y > upper)
#         3)     P(y < lower)
#         4)     E(y|lower < y < upper)
#         5)     E(y)
#         6*)    E(y|y > 0)
