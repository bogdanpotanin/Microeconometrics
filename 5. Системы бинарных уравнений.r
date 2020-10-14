# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Урок 5. Системы бинарных уравнений
# --------

# Отключим scientific notation
options(scipen = 999)

#---------------------------------------------------
# Часть 1. Система из двух бинарных уравнений
#---------------------------------------------------

#*******************************
# Пользовательский уровень
#*******************************

# Подключим дополнительные библиотеки
library("mvtnorm")                                       # симуляции из многомерного
                                                         # нормального распределения
library("numDeriv")                                      # численное дифференцирование
library("GJRM")                                          # оценивание систем
                                                         # бинарных уравнений
library("pbivnorm")                                      # двумерное нормальное распределение
 

set.seed(123)                                            # для воспроизводимости
n <- 10000                                               # число наблюдений 

# Симулируем независимые переменные
X <- rmvnorm(n,                                          # симулируем n наблюдений из многомерного
                                                         # нормального распределения
             c(0, 0, 0),                                 # с нулевым вектором математических ожиданий и
             matrix(c(1, 0.2, 0.3,                       # следующей ковариационной матрице
                      0.2, 1, -0.1,
                      0.3, -0.1, 1),
                    ncol = 3,
                    byrow = FALSE))

# Симулируем случайные ошибки
rho <- 0.3                                               # корреляция между случайными ошибками
u <- rmvnorm(n,                                          # симулируем n наблюдений из двумерного
                                                         # нормального распределения
             c(0, 0),                                    # с нулевым вектором математических ожиданий и
             matrix(c(1, rho,                            # следующей ковариационной матрице
                      rho, 1),
                    ncol = 2,
                    byrow = FALSE))

# Симулируем латентные переменные
gamma_1 <- c(1, 1.5, -2)                                 # оцениваемые регрессионные коэффициенты
gamma_2 <- c(-1, 2, -1.5)                                # для двух уравнений

z_latent_1 <- gamma_1[1] + X[, 1] * gamma_1[2] + 
                           X[, 2] * gamma_1[3] +
                           u[, 1]
z_latent_2 <- gamma_2[1] + X[, 1] * gamma_2[2] + 
                           X[, 3] * gamma_2[3] +
                           u[, 2]

# Симулируем зависимые переменные
z_1 <- as.numeric(z_latent_1 > 0)
z_2 <- as.numeric(z_latent_2 > 0)

# Сохраним все в датафрейм
d <- data.frame("z_1" = z_1, "z_2" = z_2,
                "X1" = X[, 1], "X2" = X[, 2],
                "X3" = X[, 3])

# Оцениваем параметры модели
formula_1 <- z_1 ~ X1 + X2                              # записываем формулы линейного индекса
formula_2 <- z_2 ~ X1 + X3                              # для каждого из уравнений системы
model_bp <- gjrm(formula = list(formula_1,              # задаем лист, состоящий из 
                                formula_2),             # обеих формул
                 data = d,
                 Model = "B",                           # указываем тип модели как систему
                                                        # бинанрных уравнений
                 margins = c("probit", "probit"),       # задаем маржинальные распределения
                                                        # случайных ошибок уравнений
                 BivD = "N")                            # задаем тип совместного распределения
                                                        # случайных ошибок (копулу)
                  
summary(model_bp)                                       # посмотрим на результат

# Рассмотрим корреляцию между 
# случайными ошибками уравнений
rho_est <- model_bp$theta                               # оценка корреляции между случайными ошибками
data.frame("Rho real" = rho,                            # сравним истинное значение корреляции
           "Rho estimate" = rho_est)                    # с её оценкой
cov_est <- solve(model_bp$fit$hessian)                  # оценка асимптотической ковариационной матрицы
std_rho <- sqrt(cov_est["theta.star", "theta.star"])    # оценка стандартной ошибки оценки корреляции
p_value_rho <- 2 * min(pnorm(rho_est / std_rho),        # p-value теста о равенстве корреляции между
                       1- pnorm(rho_est / std_rho))     # случайными ошибками нулю
cov_u_est <- matrix(c(1, rho_est,                       # оценка ковариационной матрицы
                      rho_est, 1),                      # совместного распределения случайных ошибок
                    ncol = 2)

# Получим оценки линейных индексов
z_li_1 <- predict(model_bp, eq = 1)                     # первое уравнение
z_li_2 <- predict(model_bp, eq = 2)                     # второе уравнение

# Оценим вероятность того, что
p_1 <- predict(model_bp,                                # значение 1                  
               type = "response",                       # будет получено
               eq = 1)                                  # в первом уравнении
p_2 <- predict(model_bp,                                # значение 1                  
               type = "response",                       # будет получено
               eq = 2)                                  # во втором уравнении

#*******************************
# Продвинутый пользовательский уровень
#*******************************

# Оценим вероятность того, что в обоих уравнениях
# будет наблюдаться успех
p_1_2 <- pbivnorm(x = cbind(z_li_1, z_li_2),            # аргументы функции распределения
                                                        # двумерного нормального распределения,
                                                        # расположенные по строкам
                  rho = rho_est)

# Оценим вероятность того, что в первом уравнении
# будет наблюдаться успех, а во втором - неудача
p_1_n2 <- pbivnorm(x = cbind(z_li_1, -z_li_2),          # аргументы функции распределения
                                                        # двумерного нормального распределения,
                                                        # расположенные по строкам
                  rho = -rho_est)

# Пользуясь формулой условной вероятности
# оценим условные вероятности того, что
p_1_cond_2 <- p_1_2 / p_2                               # в первом уравнении наблюдается успех
                                                        # при условии, что во втором также
                                                        # наблюдается успех
p_1_cond_n2 <- p_1_n2 / (1 - p_2)                       # в первом уравнении наблюдается успех
                                                        # при условии, что во втором
                                                        # наблюдается неудача

# Сравним условные и безусловные вероятности
# успеха в первом уравнении
p_df <- data.frame(p_1,                                 # безусловная вероятности
                   p_1_cond_2, p_1_cond_n2)             # условная вероятность
colnames(p_df) = c("P(z_1 = 1)", 
                   "P(z_1 = 1 | z_2 = 1)", 
                   "P(z_1 = 1 | z_2 = 0)")
head(p_df)

#*******************************
# Пользовательский уровень
#*******************************

# Рассмотрим пример на данных, отражающих готовность
# людей платить за сохранение парка
help(Kakadu)

# Загрузим данные
library("Ecdat")                                        # встроенные данные
data("Kakadu")
h <- Kakadu
h$wtp <- as.numeric(h$answer != "nn")                   # переменная принимает значение 1,
                                                        # если индивид готов заплатить за
                                                        # сохранение парка больше некоторой суммы
h$vparks <- as.numeric(h$vparks == "yes")               # переменная принимает значение 1,
                                                        # если индивид посещал парк

# Модель, описывающая готовность индивида
# заплатить более некоторой суммы
model_wtp <- glm(formula = wtp ~ age +                  # возраст
                                 sex +                  # пол (мужчина)
                                 income +               # доход в тысячах долларов
                                 moreparks +            # нужно больше парков
                                 wildlife +             # важно сохранять дикую природу
                                 aboriginal +           # важно учитывать интересы 
                                                        # коренных жителей
                                 finben,                # важно руководствоваться соображениями
                                                        # финансовой выгоды при использовании
                                                        # природных ресурсов
                 data = h,                                    
                 family = binomial(link = "probit"))          
summary(model_wtp)

# Модель, описывающая посещение индивидом парка
model_vparks <- glm(formula = vparks ~ age + sex + 
                                       schooling,               # число лет, потраченных на                                               
                    data = h,                                   # получение образования
                    family = binomial(link = "probit"))          
summary(model_vparks)

# Система бинарных уравнений
model_bp <- gjrm(list(formula(model_wtp),                       # формула второго уравнения   
                      formula(model_vparks)), 
            data=h,
            Model = "B",                                        # двумерный пробит                                      
            margins = c("probit", "probit"))                    # маржинальные распределения можно изменить
summary(model_bp)

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 1.1.     Оцените вероятность того, что индивид с
#          вашими характеристиками
#          1)     готов платить за сохранение парка
#          2)     не готов платить за сохранение парка
#          3)     готов платить за сохранение парка и
#                 не посещал этот парк
#          4)     готов платить за сохранение парка, если
#                 он не посещал этот парк
#          5*)    готов платить за сохранение парка или
#                 посещал этот парк
#          6**)   готов платить за сохранение парка, если
#                 он или готов платить за сохранение парка
#                 но не посещал парк, или не готов платить,
#                 но посещал парк
#         Примечание: используйте аргумент newdata
#                     в функции predict()
# 1.2.    На уровне значимости 5% проверьте гипотезу о независимости
#         случайных ошибок уравнений готовности платить и 
#         посещения парка
#         1)     используя стандартную ошибку оценки
#                корреляции между случайными ошибками
#         2*)    используя LR тест
# 1.3**.  Оцените предельный эффект возраста на 
#         вероятности из пункта 1.1.
#         Примечание: удобно воспользоваться
#                     численным дифференцированием
# 1.4.    На уровне значимости 5% проверьте гипотезу о том, что
#         возраст не влияет ни на вероятность проявить готовность
#         платить, ни на вероятность посещения парка
#         1)     используя LR тест
#         2)     используя тест Вальда
# 1.5*.   На уровне значимости 5% при помощи LR теста проверьте
#         возможность оценивания общей модели для мужчин и женщин
# 1.6.    Руководствуясь информационными критериями выбирете
#         1*)    оптимальные маржинальные распределения 
#                для случайных ошибок
#         2**)   оптимальную копулу и маржинальные распределения
#                для случайных ошибок
# 1.7***. Запрограммируйте процедуру оценивания системы бинарных уравнений
#         с гетероскедастичной корреляцией между случайными ошибками

#---------------------------------------------------
# Часть 2. Иерархические системы бинарных уравнений
#---------------------------------------------------

#*******************************
# Пользовательский уровень
#*******************************

# Иерархическая система подразумевает, что зависимая
# переменная одного из бинарных уравнений входит
# в качестве независимой переменной в другое
# бинарное уравнение

# Пусть независимая переменная из второго
# уравнения влияет на зависимую переменную из первого
z_latent_1 <- z_latent_1 + 2 * z_2
z_1 <- as.numeric(z_latent_1 >= 0)
d$z_1 <- z_1

# Изменим формулу для первого уравнения в соответствии
# с внесенными изменениями
formula_1 <- z_1 ~ X1 + X2 + z_2

# Оценим иерархическую систему
model_bp2 <- gjrm(formula = list(formula_1,             # задаем лист, состоящий из 
                                formula_2),             # обеих формул
                 data = d,
                 Model = "B",                           # указываем тип модели как систему
                                                        # бинанрных уравнений
                 margins = c("probit", "probit"),       # задаем маржинальные распределения
                                                        # случайных ошибок уравнений
                 BivD = "N")                            # задаем тип совместного распределения
                                                        # случайных ошибок (копулу)

summary(model_bp2)                                      # посмотрим на результат

# Применим модель к реальным данным
model_bp2 <- gjrm(list(wtp ~ age + sex + income +       # формула уравнения
                            moreparks + wildlife +      # готовности платить
                            aboriginal + finben +
                            vparks,                     # переменная на факт посещения парка     
                      vparks ~ age + sex + schooling),  # формула уравнения посещения парка
                 data=h,
                 Model = "B",                           # двумерный пробит                                      
                 margins = c("probit", "probit"))       # маржинальные распределения можно изменить
summary(model_bp2)
rho_est <- model_bp2$theta

# Создадим отдельного индивида
Boris <- data.frame(age = 30, sex = "male",
                    income = 20, moreparks = 3,
                    wildlife = 2, aboriginal = 1,
                    finben = 5,
                    vparks = 1, schooling = 8)

# Оценим линейные индексы уравнений для Бориса
z_li_1 <- predict(model_bp2, eq = 1, newdata = Boris)
z_li_2 <- predict(model_bp2, eq = 2, newdata = Boris)

# Оценим вероятности для Бориса
p_1 <- predict(model_bp2, eq = 1, type = "response",   # вероятность проявить готовность платить
               newdata = Boris)
p_2 <- predict(model_bp2, eq = 2, type = "response",   # вероятность посещать парк
               newdata = Boris)
p_1_2 <- pbivnorm(x = matrix(c(z_li_1, z_li_2),        # вероятность посещать парк и быть готовым
                             ncol = 2),                # заплатить за его сохранение
                  rho = rho_est)

#*******************************
# Продвинутый пользовательский уровень
#*******************************

# Оценим вероятность того, что Борис готов платить
# за сохранение парка при условии, что он
# посещает парк
p_1_cond_2 <- p_1_2 / p_2

# Оценим вероятность того, что Борис готов платить
# за сохранение парка при условии, что он
# не посещает парк
Boris_new <- Boris
Boris_new$vparks <- 0                                 # представим, что Борис не посещает парк
z_li_1_new <- predict(model_bp2, eq = 1,              # пересчитаем оценку линейного индекса
                      newdata = Boris_new)  
p_1_n2 <- pbivnorm(x = matrix(c(z_li_1_new,           # вероятность посещать парк и не быть 
                                -z_li_2),             # готовым заплатить за его сохранение
                              ncol = 2),               
                   rho = -rho_est)
p_1_cond_n2 <- p_1_n2 / (1 - p_2)                     # оценка искомой условной вероятности

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 2.1. Повторите задания из предыдущего раздела
#      для данной модели
# 2.2. При помощи LR теста на уровне значимости 5%
#      проверьте необходимость использования
#      иерархической модели