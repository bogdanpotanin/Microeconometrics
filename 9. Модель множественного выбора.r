# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Урок 9. Модель множественного выбора
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
library("EnvStats")                      # распределение Гумбеля 
                                         # (extreme value distribution type 1)
library("nnet")                          # мультиномиальный логит
library("mlogit")                        # мультиномиальный логит альтернативный
                                         # пакет и iia тест
library("mnlogit")                       # мультиномиальный логит еще одна 
                                         # альтернатива на которой работает IIA
                                         # из предыдущего пакета
library("mclogit")                       # мультиномиальный логит с различными
                                         # регрессорами

library("stringr")                       # работа со строками
library("tm")

# Воспроизведем процесс генерации данных,
# предполагаемый пробит моделью
# множественного выбора

# Симулируем данные
set.seed(123)                                            # для воспроизводимости
n <- 10000                                               # число наблюдений    
# Регрессоры, общие для всех альтернатив
X_c <- rmvnorm(n,                                        # симулируем n наблюдений из многомерного
                                                         # нормального распределения
             c(0, 0, 0),                                 # с нулевым вектором математических ожиданий и
             matrix(c(1, 0.2, 0.3,                       # следующей ковариационной матрице
                      0.2, 1, -0.1,
                      0.3, -0.1, 1),
                    ncol = 3,
                    byrow = FALSE))
# Регрессоры, различающиеся для альтернатив
X_d <- rmvnorm(n, rep(0, 6))
# Случайные ошибки из распределения Гумбеля
u <- cbind(revd(n), revd(n), revd(n))                                     

# Соберем в датафрейм регрессоры, общие
# для всех альтернатив
h <- data.frame("income" = X_c[, 1],                     # показатель дохода
                "health" = X_c[, 2],                     # показатель здоровья
                "age" = X_c[, 3])                        # показатель возраста
head(h, 5)
   
# Соберем в датафрейм регрессоры, различные
# для всех альтернатив
h$price.Car <- X_d[, 1] + 0.5
h$comfort.Car <- X_d[, 2] + 0.5
h$price.Taxi <- X_d[, 3] + 1
h$comfort.Taxi <- X_d[, 4] + 1
h$price.Public <- X_d[, 5]
h$comfort.Public <- X_d[, 6]
head(h, 5)
                                                   
# Создадим несколько латентных переменных, каждая
# из которых отражает предпочтения в отношении того
# или иного вида транспорта
gamma_Car <- c(0.1, 0.7, 0.6, 0.6)
gamma_Taxi <- c(0.2, 0.9, 0.8, -0.9)
gamma_Public <- c(0.3, -1, -1, 1)
gamma_d <- c(-1, 1)
z_li_Car <- gamma_Car[1] +                              # линейный индекс Машины
            h$income * gamma_Car[2] +
            h$health * gamma_Car[3] +
            h$age * gamma_Car[4] +
            h$price.Car * gamma_d[1] +
            h$comfort.Car * gamma_d[2]
z_star_Car <- z_li_Car + u[, 1]                         # латентная переменная Машины
z_li_Taxi <- gamma_Taxi[1] +                            # линейный индекс Такси
             h$income * gamma_Taxi[2] +
             h$health * gamma_Taxi[3] +
             h$age * gamma_Taxi[4] +
             h$price.Taxi * gamma_d[1] +
             h$comfort.Taxi * gamma_d[2]
z_star_Taxi <- z_li_Taxi + u[, 2]                       # латентная переменная Такси
z_li_Public <- gamma_Public[1] +                        # линейный индекс
               h$income * gamma_Public[2] +             # общественного транспорта
               h$health * gamma_Public[3] +
               h$age * gamma_Public[4] +
               h$price.Public * gamma_d[1] +
               h$comfort.Public * gamma_d[2]
z_star_Public <- z_li_Public + u[, 3]                   # латентная переменная 
                                                        # общественного транспорта

# Сформируем зависимую переменную
h$transport[(z_star_Car >= z_star_Taxi) &               # те, кто выбрал Машину
            (z_star_Car >= z_star_Public)] <- "Car"
h$transport[(z_star_Taxi >= z_star_Car) & 
            (z_star_Taxi >= z_star_Public)] <- "Taxi"   # те, кто выбрал Такси
h$transport[(z_star_Public >= z_star_Car) & 
            (z_star_Public >= z_star_Taxi)] <- "Public" # те, кто выбрал
                                                        # общественный транспорт
summary(as.factor(h$transport))

#*******************************
# Пользовательский уровень
#*******************************

# Подготовим данные
h1 <- dfidx(h,                                 # исходный датафрейм
            varying = 4:9,                     # индексы регрессоров, которые
                                               # разнятся между альтернативами
                                               # имена эти регрессоров должны 
                                               # иметь формат имя.альтернатива
            shape = "wide",                 
            choice = "transport")              # переменная, отражающая 
                                               # выбранную альтернативу

# Оценим параметры модели
model_cmlogit <- mlogit(transport ~            # формула, включающая
                        price + comfort |      # различающиеся и
                        income + health + age, # общие регрессоры
                        data = h1)             # преобразованный датафрейм
summary(model_cmlogit)   

# Сохраним оценки
coef_cmlogit <- coef(model_cmlogit,)           # все регрессионные коэффициенты
coef_Taxi <- coef_cmlogit[str_detect(          # оценки разницы в коэффициентах 
             names(coef_cmlogit), "Taxi")]     # такси и машины
coef_Public <- coef_cmlogit[str_detect(        # оценки разницы в коэффициентах
               names(coef_cmlogit), "Public")] # Общественного транспорт и машины
coef_d <- coef_cmlogit[c("price", "comfort")]  # оценки коэффициентов при
                                               # различающихся регрессорах

# Сравним истинные и предсказанные значения
data.frame("Estimate" = coef_Taxi,             # разница между коэффициентами при 
           "Real" = gamma_Taxi - gamma_Car)    # регрессорах в уравнениях Такси и Машины
data.frame("Estimate" = coef_Public,           # разница между коэффициентами при регрессорах
           "Real" = gamma_Public - gamma_Car)  # в уравнениях Такси и Машины
                                               # оценки эффектов различающихся регрессоров

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 1.1*.   Повторите приведенный пример,
#         рассмотрев собственное множество
#         из пяти альтернатив
# 1.2**.  Симулируйте случайные ошибки из совместного
#         нормального распределения и оцените
#         устойчивость результата

#---------------------------------------------------
# Часть 2. Расчет вероятностей и предельных эффектов
#---------------------------------------------------

#*******************************
# Пользовательский уровень
#*******************************

# Оценим вероятности
  # по всей выборке
probs_cmlogit <- predict(model_cmlogit,        # без аргумента newdata результат 
                         newdata = h1)         # может оказаться неправильным
head(probs_cmlogit, 5)
  # для конкретного индивида
Boris <- data.frame(income = 1,                # характеристика Бориса
                    health = 0.5,
                    age = -0.3,
                    price.Car = 0.3,
                    comfort.Car = 0.5,
                    price.Taxi = 0.7,
                    comfort.Taxi = 0.9,
                    price.Public = -0.2,
                    comfort.Public = -0.5,
                    transport = "Car")        # указываем любое из возможных
                                              # значений зависимой переменной,
                                              # что не скажется на результате
Boris1 <- dfidx(Boris, varying = 4:9,         # перекодируем Бориса в
                shape = "wide", "transport")  # нужный формат
probs_Boris <- predict(model_cmlogit,         # предскажем вероятности
                       newdata = Boris1)

# Оценим предельный эффект на вероятность
# выбора Бориса по различных переменным
  # по доходу
ME_income <- effects(model_cmlogit,
                     covariate = "income",
                     data = Boris1)
  # по цене
ME_price <- effects(model_cmlogit,
                    covariate = "price",
                    data = Boris1)

#*******************************
# Продвинутый пользовательский уровень
#*******************************

# Предельные эффекты удобно 
# рассчитывать численно
  # предельный эффект дохода
  # на вероятность каждой из альтернатив
delta <- sqrt(.Machine$double.eps)                      # приращение
Boris1_delta <- Boris1                                  # датафрейм с приращением
Boris1_delta$income <- Boris1_delta$income + delta      # добавляем приращение по
                                                        # доходу в датафрейм
probs_Boris_delta <- predict(model_cmlogit,             # считаем вероятность с
                             newdata = Boris1_delta)    # учетом приращения
ME_income_num <- (probs_Boris_delta -                   # оцениваем предельный эффект
                  probs_Boris) / delta                  # с помощью численного
                                                        # дифференцирования
data.frame("Analytical" = ME_income,                    # сравниваем численный и
           "Numeric" = ME_income_num)                   # аналитический результаты
  # предельный эффект цены поездки
  # на Машине на вероятность каждой
  # из альтернатив
delta <- sqrt(.Machine$double.eps)                      # приращение
Boris1_delta <- Boris1                                  # датафрейм с приращением
Boris1_is_car <- Boris1_delta$idx$id2 == "Car"          # строки, в которых хранится
                                                        # информация о регрессорах,
                                                        # специфических для Машины
Boris1_delta$price[Boris1_is_car] <-                    # делаем приращение цены
        Boris1_delta$price[Boris1_is_car] + delta       # поездки на Машине
probs_Boris_delta <- predict(model_cmlogit,             # считаем вероятность с
                             newdata = Boris1_delta)    # учетом приращения
ME_price_num <- (probs_Boris_delta -                    # оцениваем предельный эффект
                       probs_Boris) / delta             # с помощью численного
                                                        # дифференцирования
data.frame("Analytical" = ME_price["Car", ],            # сравниваем численный и
           "Numeric" = ME_price_num)                    # аналитический результаты

# Предельные эффекты можно посчитать и по формуле
# Рассчитеам предельный эффект на вероятность
# использования Такси
  # по доходу
ME_income_Taxi <- probs_Boris["Taxi"] * 
  (coef_cmlogit["income:Taxi"] - 
   probs_Boris["Car"] * 0 -
   probs_Boris["Taxi"] * coef_cmlogit["income:Taxi"] -
   probs_Boris["Public"] * coef_cmlogit["income:Public"])
data.frame("Analytical" = ME_income["Taxi"],               # сравним результаты
           "Numeric" = ME_income_Taxi)                     # со встроенной функцией
  # по цене поездки на общественнрм Транспорте
ME_price_Taxi_Public <- -probs_Boris["Taxi"] *
                        probs_Boris["Public"] *                  
                        coef_cmlogit["price"]
 # по цене поездки на Такси
ME_price_Taxi_Taxi <- probs_Boris["Taxi"] * 
                      (1 - probs_Boris["Taxi"]) *
                      coef_cmlogit["price"]
data.frame("Analytical" = ME_price["Taxi", "Taxi"],       # сравним результаты
           "Numeric" = ME_price_Taxi_Taxi)                # со встроенной функцией

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 2.1.    Для индивида с произвольными характеристиками
#         оцените:
#         1)     вероятность того, что он поедет на Такси
#         2)     отношение вероятностей поездки на 
#                Такси и на Машине
# 2.2.    Проверьте гипотезу о возможности исключить
#         из модели:
#         1)     Переменную на здоровье
#         2)     Все специфические для альтернатив регрессоры
#         3)     Все регрессоры, не различающиеся
#                между альтернативами
#         Подсказка: воспользуйтесь LR тестом
# 2.3.    Для индивида с произвольными характеристиками
#         оцените предельный эффект здоровья на:
#         1)     вероятность поездки на Такси
#         2*)    отношение вероятностей поездки на
#                Такси и на Машине
# 2.4.    Повторите предыдущее задание рассмотрев
#         предельный эффект комфорта
# 2.5.    Для индивида с произвольными характеристиками 
#         проверьте гипотезу о том, что он с равной
#         вероятностью:
#         1**)   Поедет на Такси или на Машине
#         2**)   Поедет на любом из видов транспорта
# 2.6.    Постройте 95% доверительный интервал для
#         предельного эффекта:
#         1**)  цены поездки на Такси на вероятность
#               поездки на Машине
#         2**)  цены поездки на Такси на вероятность
#               поездки на Такси
#         3**)  здоровья индивид на вероятность
#               поездки на Такси

#---------------------------------------------------
# Часть 3. Тестирование гипотезы о независимости
#          от посторонних альтернатив IIA
#---------------------------------------------------

mlogit::hmftest(model_cmlogit, z = c("Car", "Taxi"))    # независимость выбора между
                                                        # Машиной и Такси от возможности
                                                        # использовать Общественный транспорт
mlogit::hmftest(model_cmlogit, z = c("Car", "Public"))  # независимость выбора между
                                                        # Машиной и Общественным транспортом 
                                                        # от возможности использовать Такси
mlogit::hmftest(model_cmlogit, z = c("Taxi", "Public")) # независимость выбора между
                                                        # Такси и Общественным транспортом
                                                        # от возможности использовать Машину

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 3.1*.   Повторите пример с пятью альтернативами и
#         проверьте возможность удалить любые две 
#         из них

#---------------------------------------------------
# Часть 4. Множественная пробит модель
#---------------------------------------------------

# Очень долго оцениваем параметры модели
model_cmprobit <- mlogit(transport ~                    # синтаксис такой же, как    
                         price + comfort |              # для логит модели
                         income + health + age, 
                         data = h1, 
                         probit = TRUE)                 # добавляем этот аргумент            
summary(model_cmprobit)   

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 4.1***. Дождитесь окончания расчетов пробит модели
#         множественного выбора