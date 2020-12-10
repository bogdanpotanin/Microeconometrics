# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Урок 12. Модель Ньюи и модель Галланта и Нички
# --------

# Отключим scientific notation
options(scipen = 999)

# Подключим дополнительные библиотеки
library("mvtnorm")                       # симуляции из многомерного
                                         # нормального распределения

library("numDeriv")                      # численное дифференцирование

library("stringr")                       # работа со строками
library("tm")

library("sampleSelection")               # метод Хекмана

library("hpa")                           # метод Галланта и Нички

#---------------------------------------------------
# Часть 1. Оценивание параметров
#---------------------------------------------------
set.seed(123)
# Зададим количество наблюдений
n <- 5000
# Истинные значения регрессионных коэффициентов
beta <- c(1, 2, 3)
gamma <- c(0.5, 1, 2)
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
# Симулируем случайные ошибки из смеси
# распределений Cтьюдента
df <- 5                                        # зададим степени свободы
e_1 <- rt(n, df)                               # симулируем выборки из трех
e_2 <- rt(n, df)                               # независимых распределений
e_3 <- rt(n, df)                               # Стьюдента с df степенями свободы
e_b <- rbinom(n, 1, 0.5)                       # симулируем выборку из
                                               # распределения Бернулли
e_d <- 3                                       # разница между математическими
                                               # ожиданиями компонент смеси
epsilon <- matrix(0, ncol = 2, nrow = n)       # итоговые случайные ошибки
epsilon[, 1] <- (e_1 + e_2 + e_d) * e_b + 
                (e_3 - e_d) * (1 - e_b)
epsilon[, 2] <- (e_1 + e_3 + e_d) * e_b + 
                (e_2 - e_d) * (1 - e_b)
rho <- cor(epsilon)[1, 2]                      # настоящая около 0.855, с малой
                                               # погрешностью будем предполагать,
                                               # что rho совпадает с данной оценкой
var(epsilon[,2])                               # настоящая около 11.5
plot(density(epsilon[, 1]))                    # визуализируем распределение

#*******************************
# Пользовательский уровень
#*******************************

set.seed(123)

# Латентная зависимая переменная
  # уравнения отбора 
z_star <- gamma[1] + gamma[2] * X[, 1] + 
          gamma[3] * X[, 2] + epsilon[, 1]
  # основного уравнения
y_star <- beta[1] + beta[2] * X[, 2] + 
          beta[3] * X[, 3] + epsilon[, 2]

# Наблюдаемая зависимая переменная
  # основного уравнения
z <- ifelse(z_star >= 0, 1, 0)
  # уравнения отбора 
y <- ifelse(z == 1, y_star, NA)

# Сюжет:
# Представим, что вы изучаете, как различные 
# факторы влияют на затраты индивидов на содержание 
# кота. При этом вы наблюдаете объем затрат лишь для
# индивидов, у которых есть кот.
h <- data.frame("cost" = y,                          # затраты на кота
                "cat" = z,                           # факт наличия кота
                "age" = X[, 1],                      # возраст
                "hobby" = X[, 2],                    # увлечение котами
                "income" = X[, 3])                   # доход

# Оценим параметры модели при помощи
  # МНК
model_ls <- lm(cost ~ hobby + income,  
               data = h)
summary(model_ls)                                    # результат оценивания
coef_ls <- coef(model_ls)                            # сохраним оценки коэффициентов
sigma(model_ls)                                      # оценка стандартного отклонения
  # метода Хекмана, основанный на ММП
model_mle <- selection(                              
  selection = cat ~ age + hobby,                     # уравнение отбора
  outcome = cost ~ hobby + income,                   # основное уравнение
  data = h,                                          # данные
  method = "ml")                                     # метод расчета ММП
summary(model_mle)                                   # результат оценивания
coef_mle <- coef(model_mle, part = "outcome")        # сохраним оценки коэффициентов
rho_mle <- model_mle$estimate["rho"]                 # оценка корреляции между
                                                     # случайными ошибками
sigma_mle <- model_mle$estimate["sigma"]             # стандартное отклонение
                                                     # случайной ошибки
  # метода Хекмана, основанного на
  # двухшаговой процедуре
model_2st <- selection(                              
  selection = cat ~ age + hobby,                     
  outcome = cost ~ hobby + income,                   
  data = h,                                          
  method = "2step")                                  # метод расчета двухшаговая процедура
summary(model_2st)                                   # результат оценивания
coef_2st <- coef(model_2st, part = "outcome")        # сохраним оценки коэффициентов
coef_2st <- coef_2st[-length(coef_2st)]              # удалим лишний коэффициент
rho_2st <- model_2st$rho                             # оценка корреляции между
                                                     # случайными ошибками

# Метод Галланта и Нички
model_gn <- hpaSelection(                              
  selection = cat ~ age + hobby,                     # формула                   
  outcome = cost ~ hobby + income,                         
  data = h,                                          # данные
  is_Newey_loocv = TRUE,                             # кросс-валидация при помощи
                                                     # критерия leave-one-out
  pol_elements = 5,                                  # число степеней полинома
                                                     # в используемом на первом шаге
                                                     # методе Ньюи
  z_K = 2,                                           # степени полинома для аппроксимации
  y_K = 2)                                           # случайных ошибок
summary(model_gn)
coef_gn <- model_gn$y_coef                           # оценки коэффициентов
                                                     # основного уравнения
coef_gn <- c( 
  model_gn$re_moments$selection_exp, coef_gn)        # оценка константы как оценка
names(coef_gn)[1] <- "(Intercept)"                   # безусловного математического
                                                     # ожидания случайной ошибки
rho_gn <- model_gn$re_moments$rho                    # оценка корреляции между
                                                     # случайными ошибками
sigma_gn <-model_gn$re_moments$selection_var         # стандартное отклонение
                                                     # случайной ошибки
# Визуализируем оценку плотности случайных ошибок
plot(model_gn, is_outcome = TRUE)                    # основное уравнение
plot(model_gn, is_outcome = FALSE)                   # уравнение отбора

# Метод Ньюи
# Автоматически применяется на первом шаге перед
# оцениванием метода Галланта и Нички
model_nw <- model_gn$Newey                               
coef_nw <- c(model_nw$constant_biased,               # оценки коэффициентов
             model_nw$y_coef)
names(coef_nw) <- names(coef_gn)

# Сравним оценки и истинные значения
  # регрессионные коэффициенты
data.frame("Real" = beta,                            # истинные значения
           "Least Squares" = coef_ls,                # МНК оценки
           "Heckman MLE" = coef_mle,                 # оценки ММП Хекмана
           "Heckman 2step" = coef_2st,               # оценки двухшагового Хекмана
           "Gallant and Nychka" = coef_gn,           # оценки Галланта и Нички 
           "Newey" = coef_nw)
  # корреляция случайных ошибок
data.frame("Real" = rho,                             # истиннoе значение
           "Least Squares" = 0,                      # МНК оценки
           "Heckman MLE" = rho_mle,                  # оценка ММП Хекмана
           "Heckman 2step" = rho_2st,                # оценка двухшагового Хекмана
           "Galland and Nychka" = rho_gn)            # оценка метода Галланти и Нички

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 1.1*.   Придумайте и реализуйте собственный
#         пример на симулированных данных


#---------------------------------------------------
# Часть 2. Расчет предсказаний и предельных эффектов
#---------------------------------------------------

#*******************************
# Пользовательский уровень
#*******************************

# Создадим индивида
Boris <- data.frame("cost" = 1,
                    "cat" = 1,
                    "age" = 0.1,
                    "hobby" = 0.3,
                    "income" = 0.5)

# Рассчитаем оценку безусловного математического 
# ожидания зависимой переменной основного уравнения,
# то есть E(y*)
  # настоящее значение
cost_star <- beta[1] + beta[2] * Boris$hobby + 
                       beta[3] * Boris$income
  # метод Хекмана
cost_star_mle <- predict(model_mle, 
                         newdata = Boris, 
                         part = "outcome",                   # для основного уравнения
                                                             # строится предсказание
                         type = "unconditional")             # безусловные предсказания  
  # метод Галланта и Нички
cost_star_gn <- predict(model_gn,
                        newdata = Boris,
                        is_cond = FALSE,                     # безулоснове предсказание
                        is_outcome = TRUE)                   # для основного уравнения
  # метод Ньюи
cost_star_nw <- predict(model_gn,
                        newdata = Boris,
                        is_cond = FALSE,                     # безусловное предсказание
                        is_outcome = TRUE,                   # для основного уравнения
                        method = "Newey")$y +                # метод Ньюи  
                model_gn$Newey$constant_biased               # и в конце прибавляем
                                                             # смещенную оценку константы
# Сравнение
data.frame("Real" = cost_star,
           "Heckman" = cost_star_mle,
           "Gallant and Nychka" = as.numeric(cost_star_gn),
           "Newey" = cost_star_nw)

# Рассчитаем оценку условного математического 
# ожидания зависимой переменной основного уравнения,
# то есть E(y*|z)
  # метод Хекмана
cost_cond_mle <- predict(model_mle, 
                         newdata = Boris, 
                         part = "outcome",               # для основного уравнения
                         type = "conditional")           # условные предсказания 
cost_cond_mle[1]                                         # E(y*|z = 0)
cost_cond_mle[2]                                         # E(y*|z = 1)
  # метод Галланта и Нички
cost_cond_gn <- predict(model_gn,
                        newdata = Boris,
                        is_cond = TRUE,                  # условное предсказание
                        is_outcome = TRUE)               # для основного уравнения
cost_cond_gn$y_0[1]                                      # E(y*|z = 0)
cost_cond_gn$y_1[1]                                      # E(y*|z = 1)
  # метод Ньюи
cost_cond_nw <- predict(model_gn,
                        newdata = Boris,
                        is_cond = TRUE,                   # безусловное предсказание
                        is_outcome = TRUE,                # для основного уравнения
                        method = "Newey")                 # метод Ньюи  
# Сравнение
  # E(y*|z = 0)
data.frame("Heckman" = cost_cond_mle[1] ,
           "Gallant and Nychka" = cost_cond_gn$y_0[1])
  # E(y*|z = 1)
data.frame("Heckman" = cost_cond_mle[2],
           "Gallant and Nychka" = cost_cond_gn$y_1[1],
           "Newey" = cost_cond_nw)

# Оценим P(z = 1)
  # метод Хекмана
cat_prob_mle <- predict(model_mle, 
                        newdata = Boris, 
                        part = "selection",              # для уравнения отбора
                        type = "response")               # предсказываем вероятность
  # метод Галланта и Нички
cat_prob_gn <- predict(model_gn,
                       newdata = Boris,
                       is_cond = FALSE,                 # безусловное предсказание
                       is_outcome = FALSE)              # для уравнения занятости

# Оценим линейный индекс
  # Методом Хекмана
cat_li <- predict(model_mle, 
                  newdata = Boris, 
                  part = "selection",                    # для уравнения отбора
                  type = "link")                         # предсказываем линейный индекс
lambda_est_1 <- dnorm(cat_li) / pnorm(cat_li)            # оценка отношения Миллса

#*******************************
# Продвинутый пользовательский уровень
#*******************************

# Рассчитаем предельный эффект увлечения котами
# на затраты на содержание кота на E(y*|z=1)

coef_s_est <- coef(model_mle)[1:3]
  # аналитически 
    # методом Хекмана
hobby_ME_mle <- coef_mle["hobby"] - rho_mle * sigma_mle *
                                    (cat_li * lambda_est_1 +
                                    lambda_est_1 ^ 2) *
                                    coef_s_est["hobby"]
    # методом Галланта и Нички
      # создадим индивида с приращением
eps <- 1e-8
Boris_eps <- Boris
Boris_eps$hobby <- Boris$hobby + eps
      # осуществим расчет
cost_cond_gn_eps <- predict(model_gn,
                            newdata = Boris_eps,
                            is_cond = TRUE,              # условное предсказание
                            is_outcome = TRUE)           # для основного уравнения
hobby_ME_gn <- (cost_cond_gn_eps$y_1 - 
                cost_cond_gn$y_1) / eps
  # методом Ньюи
cost_cond_nw_eps <- predict(model_gn,
                            newdata = Boris_eps,
                            is_cond = TRUE,                   # безусловное предсказание
                            is_outcome = TRUE,                # для основного уравнения
                            method = "Newey")                 # метод Ньюи  
hobby_ME_nw <- (cost_cond_nw_eps$y_1 - 
                  cost_cond_nw$y_1) / eps

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 1.1.    Для индивида с произвольными характеристиками
#         оцените с помощью метода Ньюи и метода
#         Галланта и Нички предельный эффект на:
#         1)     E(y*)
#         2)     P(z = 1)
#         3)     E(y|z = 1)
#         4)     E(y|z = 0)
#         5*)    P(z = 1|y = 0.5)
# 1.2.    Повторите предыдущее задание добавив
#         в модель:
#         1*)    доход в квадрате
#         2*)    взаимодействие между возрастом
#                и доходом

#---------------------------------------------------
# Часть 3. Метод Ньюи
#---------------------------------------------------

#*******************************
# Продвинутый пользовательский уровень
#*******************************

# Напишем вспомогательную функции для подсчета MSE 
# при leave-one-out кросс валидации
loocv <- function(fit)
{
  h <- lm.influence(fit)$h
  loocv_value <- mean((residuals(fit)/(1 - h)) ^ 2)
  
  return(loocv_value)
}

# Напишем функцию для расчета отношения Миллса
# стандартного нормального распределения
mills <- function(x)
{
  value <- dnorm(x) / pnorm(x)
  
  return(value)
}

# Реализуем самостоятельно метод Ньюи

# На первом шаге
model_1 <- hpaBinary(cat ~ age + hobby,                 # оценим уравнение отбора
                     data = h,                          # методом Галланта и Нички
                     K = 3)
summary(model_1)
plot(model_1)                                           # визуализируем результат
cat_li <- model_1$z_latent                              # достанем линейный индекс
h$lambda <- mills((cat_li + model_1$errors_exp) /       # оценим лямбду
                  sqrt(model_1$errors_var))          

# На втором шаге оценим
  # модель с первой степенью
model_2_1 <- lm(cost ~ hobby + income +                 # используем МНК включив
                lambda,         
              data = h)                                 # оценку лямбды как регрессор                           
coef_2_1 <- coef(model_2_1)[1:3]                        # достаем оценки коэффициентов
loocv_2_1 <- loocv(model_2_1)                           # сохраним leave-one-out RMSE
  # модель со второй степенью
model_2_2 <- lm(cost ~ hobby + income +                 # используем МНК включив
                  lambda + I(lambda ^ 2),         
                data = h)                               # оценку лямбды как регрессор                           
coef_2_2 <- coef(model_2_2)[1:3]                        # достаем оценки коэффициентов
loocv_2_2 <- loocv(model_2_2)                           # сохраним leave-one-out RMSE
  # модель с третьей степенью
model_2_3 <- lm(cost ~ hobby + income +                 # используем МНК включив
                  lambda + 
                  I(lambda ^ 2) + I(lambda ^ 3),         
                data = h)                               # оценку лямбды как регрессор                           
coef_2_3 <- coef(model_2_3)[1:3]                        # достаем оценки коэффициентов
loocv_2_3 <- loocv(model_2_3)                           # сохраним leave-one-out RMSE

# Сопоставим результаты
c(loocv_2_1, loocv_2_2, loocv_2_3)                      # по leave-one-out RMSE
data.frame("Real" = beta,
           "Model 1" = coef_2_1,
           "Model 2" = coef_2_2,
           "Model 3" = coef_2_3)

# Бутстрап для получения ковариационной
# матрицы регрессионных коэффициентов
boot_iter <- 20                                           # число бутстрап итераций, на практике
                                                          # используется хотя бы 100
coef_boot <- matrix(NA, nrow = boot_iter, ncol = 3)       # сохраняем коэффициенты
                                                          # каждую бутстрап итерацию
colnames(coef_boot) <- names(coef_2_1)
for(i in 1:boot_iter)
{
  boot_ind <- sample(1:n, n, replace = TRUE)              # генерируем новую выборку
  h_boot <- h[boot_ind, ]                                 # с возвращением
  # На первом шаге
  model_boot_1 <- hpaBinary(cat ~ age + hobby,            # оценим уравнение отбора
                            data = h_boot,                # методом Галланта и Нички
                            K = 3,
                            x0 = model_1$x1)              # для ускорения каждый раз стартуем
                                                          # с хорошей точки
  cat_li_boot <- model_boot_1$z_latent                    # достанем линейный индекс
  h_boot$lambda <- mills((cat_li_boot + 
                          model_boot_1$errors_exp) /      # оценим лямбду
                          sqrt(model_boot_1$errors_var))          
  
  # На втором шаге оценим
  # модель с первой степенью
  model_2_1_boot <- lm(cost ~ hobby + income +            # используем МНК включив
                       lambda,         
                       data = h_boot)                     # оценку лямбды как регрессор                           
  coef_2_1_boot <- coef(model_2_1_boot)[1:3]              # достаем оценки коэффициентов
  loocv_2_1_boot <- loocv(model_2_1_boot)                 # сохраним leave-one-out RMSE
  # модель со второй степенью
  model_2_2_boot <- lm(cost ~ hobby + income +            # используем МНК включив
                       lambda + I(lambda ^ 2),         
                       data = h_boot)                     # оценку лямбды как регрессор                           
  coef_2_2_boot <- coef(model_2_2_boot)[1:3]              # достаем оценки коэффициентов
  loocv_2_2_boot <- loocv(model_2_2_boot)                 # сохраним leave-one-out RMSE
  # модель с третьей степенью
  model_2_3_boot <- lm(cost ~ hobby + income +            # используем МНК включив
                       lambda + 
                       I(lambda ^ 2) + I(lambda ^ 3),         
                       data = h_boot)                     # оценку лямбды как регрессор                           
  coef_2_3_boot <- coef(model_2_3_boot)[1:3]              # достаем оценки коэффициентов
  loocv_2_3_boot <- loocv(model_2_3_boot)                 # сохраним leave-one-out RMSE
  
  coef_all <- cbind(coef_2_1_boot, coef_2_2_boot, 
                    coef_2_3_boot)                        # сохраняем все коэффициенты
  loocv_min <- which.min(c(loocv_2_1_boot,                # находим модель с меньшим RMSE
                           loocv_2_2_boot, 
                           loocv_2_3_boot))
  
  coef_boot[i, ] <- coef_all[, loocv_min]                 # выбираем оценки модели с
                                                          # наименьшим RMSE
}
# Получаем и используем состоятельную
# оценку ковариационной матрицы
library("lmtest")
cov_boot <- cov(coef_boot)
c(quantile(coef_boot[, "hobby"], 0.05),                   # 90% бутстрапированный ДИ
  quantile(coef_boot[, "hobby"], 0.95))

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 1.1*.   Проверьте значимость регрессоров ориентируясь
#         на асимптотическую нормальность оценок 
#         метода Ньюи
# 1.2**.  Проверьте гипотезу о возможности оценивания
#         МНК модели вместо метода Ньюи