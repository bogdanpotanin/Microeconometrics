# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Урок 4. Тестирование гипотез о распределении 
#         случайных ошибок моделей бинарного выбора
# --------

# Отключим scientific notation
options(scipen = 999)

#---------------------------------------------------
# Часть 1. Учет и тестирование гетероскедастичности
#          в моделях бинарного выбора
#---------------------------------------------------

#*******************************
# Пользовательский уровень
#*******************************

library("titanic")                                               # пакет, содержащий датафрейм с
                                                                 # информацией о пассажирах Титаника

library("glmx")                                                  # пакет, позволяющий оценивать пробит
                                                                 # модель с гетероскдестичной 
                                                                 # случайной ошибкой

library("lmtest")                                                # дополнительные тесты

library("numDeriv")                                              # численное дифференцирование

# Достаем данные о выживших на титанике 
# из встроенного набора данных
data(titanic_train)                                              # загружаем встроенный датафрейм
h <- na.omit(titanic_train)                                      # сохраняем его
help(titanic_train)                                              # и читаем сопутствующую документацию

# Преобразуем некотонрые переменные и во избежание
# проблем при дальнейшем анализе преобразуем все
# переменные к численному типу
h$sibl <- h$SibSp + h$Parch                                      # колиечство близких на борту
h$class_1 <- as.numeric(h$Pclass == 1)                           # пассажир плывет первым классом
h$class_2 <- as.numeric(h$Pclass == 2)                           # пассажир плывет вторым классом
h$male <- as.numeric(h$Sex == "male")                            # пассажир мужчина

# Построим пробит модель для того чтобы оценить, как
# различные характеристики индивидов влияли на вероятность
# того, что они выживут на Титанике
model_probit <- glm(formula = Survived ~ Age + I(Age ^ 2) +      # указываем формулу без константы, поскольку
                                         male + sibl +           # она учитывается автоматически
                                         class_1 + class_2,      
                    data = h,                                    # датафрейм, из которого берутся зависимая
                                                                 # и независимые переменные
                    family = binomial(link = "probit"))          # тип оцениваемой бинарной регрессии: в данном
summary(model_probit)

# Теперь учтем возможное различие в дисперсия случайных
# ошибок, обусловленное полом и возрастом
model_hetprobit <- hetglm(formula = Survived ~ Age + I(Age ^ 2) +   # основное уравнение
                                               male + sibl +                   
                                               class_1 + class_2 |  # после | следует  
                                               Age + male,          # уравнение линейного индекса 
                                                                    # дисперсии случайной ошибки      
                          data = h,                                 # датафрейм, из которого берутся зависимая
                                                                    # и независимые переменные
                          family = binomial(link = "probit"),       # распределение случайных ошибок
                          link.scale = "log")                       # тип функции, которая берется от линейного 
                                                                    # индекса уравнения дисперсии случайных ошибок
summary(model_hetprobit)
# Обратите внимание, что в формуле после | стоят регрессоры, 
# влияющие на дисперсию случайной ошибки. Обозначим линейный
# индекс уравнения дисперсии как tau * W, где:
# W - матрица независимых переменных, влияющих на sigma
# tau - коэффициенты при W
# В функции hetglm() link scale указывает, 
# в каком виде представлена ошибка. 
# Имеются следующие варианты:
# 1. identity  ---  sigma        =  tau * W
# 2. log       ---  log(sigma)   =  tau * W  =>  sigma = exp(tau * W)
# 3. sqrt      ---  sqrt(sigma)  =  tau * W  =>  sigma = (tau * W) ^ 2
# Примечание: для краткости индекс наблюдений пропущен

# Осуществим тест на гомоскедастичность:
# H0: tau = 0
lrtest(model_hetprobit, model_probit) 

# Сравним вероятности, предсказанные двумя моделями
p_probit <- predict(model_probit, type = "response")                # оценки вероятности при допущении о
                                                                    # гомоскедастичности случайных ошибок
p_hetprobit <- predict(model_hetprobit, type = "response")          # оценки вероятности по модели с
                                                                    # гетероскедастичной случайной ошибкой
plot(p_probit, p_hetprobit,                                         # визуально сравним оценки
     main = "Probabilities estimates",
     xlab = "Probit",
     ylab = "Heteroscedastic probit")

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 1.1.    Оцените модель, в которой на дисперсию
#         случайной ошибки влияет класс, которым
#         плывет пассажир. Проверьте наличие
#         гетероскедастичности.
# 1.2.    Сравните модели с учетом гетероскедастичности
#         случайной ошибки и без нее с использованием
#         информационный критериев AIC и BIC
#         Подсказка: используйте функции AIC() и BIC
#         в отношении возвращаемых функциями glm()
#         и hetglm() объектов


#*******************************
# Академический уровень
#*******************************

# Достаем оценки коэффициентов
gamma_est <- matrix(model_hetprobit$coefficients$mean,              # оценки gamma
                    ncol = 1)
rownames(gamma_est) <- names(model_hetprobit$coefficients$mean)
tau_est <- matrix(model_hetprobit$coefficients$scale,
                  ncol = 1)                                         # оценки tau
rownames(tau_est) <- names(model_hetprobit$coefficients$scale)
# Достанем матрицы независимых переменных
X <- model.frame(model_hetprobit)
X$'(Intercept)' <- 1
X_names <- colnames(X)

# Получим оценки линейного индекса 
z_li_mean <- as.matrix(X[, rownames(gamma_est)]) %*% gamma_est
z_li_scale <- as.matrix(X[, rownames(tau_est)]) %*% tau_est

# Оценим вероятности выжить
p_hetprobit <- pnorm(z_li_mean,
                     sd = exp(z_li_scale))

# Оценим предельный эффект возраста
# на вероятность выжить
age_ME <- dnorm(z_li_mean,
                sd = exp(z_li_scale)) *
          ((gamma_est["Age", ] + 
            2 * gamma_est["I(Age^2)", ] * h$Age) -
          tau_est["Age", ] * z_li_mean) / 
          exp(z_li_scale)

# Проверим гипотезу о гетероскедастичности
# с помощью LM теста

# Запишем функцию правдоподобия для данной модели
HetprobitLnL <- function(x, z, X, W, scale_fn,        # функция правдоподобия
                         is_aggregate = TRUE)                   
{
  m_X <- ncol(X)
  m_W <- ncol(W)
  
  gamma <- matrix(x[1:m_X], ncol = 1)                 # вектор gamma коэффициентов и
  tau <- matrix(x[(m_X + 1):(m_X + m_W)], ncol = 1)   # вектор дополнительных параметров  
                                                      # переводим в матрицу с одним столбцом

  z_li_mean <- X %*% gamma                            # оценка линейного индекса
  z_li_scale <- W %*% tau                             # латентной переменной
  z_li_scale_fn <- scale_fn(z_li_scale)
  
  n_obs <- nrow(X)                                       # количество наблюдений
  
  L_vec <- matrix(NA, nrow = n_obs,                      # вектор столбец вкладов наблюдений
                  ncol = 1)                              # в функцию правдоподобия
  
  is_z_0 <- z == 0                                       # вектор условий z = 0
  is_z_1 <- z == 1                                       # вектор условий z = 1

  L_vec[is_z_1] <- pnorm(z_li_mean[is_z_1], 
                         sd = z_li_scale_fn[is_z_1])     # вклад наблюдений для которых zi = 1
  L_vec[is_z_0] <- 1 - pnorm(z_li_mean[is_z_0],
                             sd = z_li_scale_fn[is_z_0]) # вклад наблюдений для которых zi = 0
  
  lnL_vec <- log(L_vec)                                  # логарифмы вкладов
  
  if(!is_aggregate)                                      # возвращаем вклады
  {                                                      # при необходимости
    return(lnL_vec)
  }
  
  lnL <- sum(lnL_vec)                                    # логарифм функции правдоподобия
  
  return(lnL)
}

# Достанем данные
df_hetprobit <- model.frame(model_hetprobit, )
X_mat <- as.matrix(df_hetprobit)
X_mat[, 1] <- 1
W_mat <- as.matrix(X_mat[, c("Age", "male")])

# Достанем оценки обычной модели
x_est_R <- c(model_probit$coefficients, 
             rep(0, ncol(W_mat)))

# Применим функцию
lnL_R <- HetprobitLnL(x_est_R, df_hetprobit$Survived,
                      X_mat, W_mat, exp)                 # считаем логарифм функции правоподобия
# при ограничениях, совпадающую с логарифмом
# функции правдоподобия обычной пробит модели
lnL_R_grad <- grad(func = HetprobitLnL,                  # считаем градиент данной функции
                   x = x_est_R,                          # численным методом
                   z = df_hetprobit$Survived, 
                   X = X_mat, W = W_mat,
                   scale_fn = exp)                       # замените exp на function(x)
                                                         #                 {
                                                         #                   return(abs(x + 1)})
                                                         #                 }
                                                         # и убедитесь, что результат не изменится
lnL_R_grad <- matrix(lnL_R_grad, ncol = 1)               # градиент как матрица с одним столбцом
lnL_R_Jac <- jacobian(func = HetprobitLnL,               # оцениваем асимптотическую ковариационную
                   x = x_est_R,                          # матрицу при помощи Якобиана, расчитанного
                   z = df_hetprobit$Survived,            # численным методом, поскольку численно
                   X = X_mat, W = W_mat,                 # рассчитать Гессиан достаточно точным 
                   scale_fn = exp,                       # образом не получается
                   is_aggregate = FALSE)
as_cov_est <- solve(t(lnL_R_Jac) %*% lnL_R_Jac)          # cчитаем оценку асимптотической ковариационной
                                                         # матрицы с помощью Якобианов, поскольку 
                                                         # численным методом Гессиан считается
                                                         # очень плохо
# Реализуем тест
LM_value <- t(lnL_R_grad) %*%                            # считаем статистику теста
  as_cov_est %*%                                         # множителей Лагранжа
  lnL_R_grad
p_value <- 1 - pchisq(LM_value, df = 2)                  # рассчитываем p-value теста

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 1.1.    Убедитесь, что статистика теста не изменится,
#         если в качестве функции, определяющей дисперсию
#         случайной ошибки, будет выступать квадрат суммы
#         линейного индекса (уравнения дисперсии) и единицы.
#         Привидите собственный пример функции, при которой
#         статистика теста останется прежней.
# 1.2*.   Реализуйте LM тест на гетероскедастичность для
#         логит модели

#---------------------------------------------------
# Часть 2. Выбор между моделями бинарного выбора
#          с различным распределением случайных
#          ошибок на основе информационных
#          критериев
#---------------------------------------------------

#*******************************
# Пользовательский уровень
#*******************************

# Оценим несколько моделей с альтернативным
# распределением случайных ошибок
  
# Логистическое распределение
model_logit <- glm(formula = Survived ~ Age + I(Age ^ 2) +       # указываем формулу без константы, поскольку
                             male + sibl +                       # она учитывается автоматически
                             class_1 + class_2,      
                   data = h,                                     # датафрейм, из которого берутся зависимая
                                                                 # и независимые переменные
                   family = binomial(link = "logit"))            # тип оцениваемой бинарной регрессии: в данном
summary(model_logit)

# Распределение Коши
model_cauchit <- glm(formula = Survived ~ Age + I(Age ^ 2) +     # указываем формулу без константы, поскольку
                               male + sibl +                     # она учитывается автоматически
                               class_1 + class_2,      
                     data = h,                                   # датафрейм, из которого берутся зависимая
                                                                 # и независимые переменные
                   family = binomial(link = "cauchit"))          # тип оцениваемой бинарной регрессии: в данном
summary(model_cauchit)

# Для самых любопытных
# 
# library("hpa")
#
# Распределение Галланта и Нички с
# четвертым порядком полинома
# model_GN <- hpaBinary(formula = Survived ~ class_1 + class_2 +
#                                 male + sibl +                     
#                                 Age + I(Age ^ 2),      
#                       data = h, 
#                       K = 4, 
#                       cov_type = "sandwichFD")                                   
# summary(model_GN)
# AIC(model_GN)
# plot(model_GN)

# Сравним модели по информационным 
# критериям AIC и BIC
  # AIC
AIC_probit <- AIC(model_probit)
AIC_logit <- AIC(model_logit)
AIC_cauchit <- AIC(model_cauchit)
  # BIC
BIC_probit <- BIC(model_probit)
BIC_logit <- BIC(model_logit)
BIC_cauchit <- BIC(model_cauchit)
  # сравнение
IC_df <- data.frame("AIC" = c(AIC_probit, 
                              AIC_logit, 
                              AIC_cauchit),
                    "BIC" = c(BIC_probit, 
                              BIC_logit, 
                              BIC_cauchit))
rownames(IC_df) <- c("Normal", 
                     "Logistic", 
                     "Cauchit")
print(IC_df)

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 1.1.    Удалите переменные на возраст и проверьте, какая
#         из моделей окажется предпочтительней согласно
#         критериям AIC и BIC

#---------------------------------------------------
# Часть 3. Тестирование гипотезы о
#          нормальном распределении
#---------------------------------------------------

#*******************************
# Академический уровень
#*******************************

# Запишем функцию правдоподобия
# для модели со случайно ошибкой
# из распределения Пирсона
ProbitLnLExtended <- function(x, z, X,                   # функция правдоподобия
                              is_aggregate = TRUE)                   
{
  gamma <- matrix(x[-c(1, 2)], ncol = 1)                 # вектор gamma коэффициентов и
  t <- matrix(x[c(1, 2)], ncol = 1)                      # вектор дополнительных параметров  
                                                         # переводим в матрицу с одним столбцом
  z_li <- X %*% gamma                                    # оценка линейного индекса
  z_est <- z_li + t[1] * z_li ^ 2 +                      # оценка математического ожидания 
                  t[2] * z_li ^ 3                        # латентной переменной
  
  n_obs <- nrow(X)                                       # количество наблюдений
  
  L_vec <- matrix(NA, nrow = n_obs,                      # вектор столбец вкладов наблюдений
                  ncol = 1)                              # в функцию правдоподобия
  
  is_z_0 <- z == 0                                       # вектор условий z = 0
  is_z_1 <- z == 1                                       # вектор условий z = 1
  
  L_vec[is_z_1] <- pnorm(z_est[is_z_1])                  # вклад наблюдений для которых zi = 1
  L_vec[is_z_0] <- 1 - pnorm(z_est[is_z_0])              # вклад наблюдений для которых zi = 0
  
  lnL_vec <- log(L_vec)                                  # логарифмы вкладов
  
  if(!is_aggregate)                                      # возвращаем вклады
  {                                                      # при необходимости
    return(lnL_vec)
  }
  
  lnL <- sum(lnL_vec)                                    # логарифм функции правдоподобия
  
  return(lnL)
}
# Воспользуемся созданной функцией
  # Оценки модели при справедливом ограничении,
  # накладываемом нулевой гипотезой
gamma_est <- coef(model_probit)                          # достаем оценки из обычной пробит
gamma_R <- c(0, 0, gamma_est)                            # модели и приравниваем значения
names(gamma_R)[c(1, 2)] <- c("t1", "t2")                 # дополнительных параметров к значениям,
                                                         # предполагаемым нулевой гипотезой
  # Создадим матрицу регрессоров
X_mat <- as.matrix(model.frame(model_probit))            # достаем датафрейм с регрессорами и
X_mat[, 1] <- 1                                          # первращаем его в матрицу, а также
colnames(X_mat)[1] <- "Intercept"                        # заменяем зависимую переменную на константу
  # Применим функцию
lnL_R <- ProbitLnLExtended(gamma_R, h$Survived, X_mat)   # считаем логарифм функции правоподобия
                                                         # при ограничениях, совпадающую с логарифмом
                                                         # функции правдоподобия обычной пробит модели
lnL_R_grad <- grad(func = ProbitLnLExtended,             # считаем градиент данной функции
                   x = gamma_R,                          # численным методом
                   z = h$Survived, X = X_mat)
lnL_R_grad <- matrix(lnL_R_grad, ncol = 1)               # градиент как матрица с одним столбцом
lnL_R_Jac <- jacobian(func = ProbitLnLExtended,          # считаем Якобин данной функции
                      x = gamma_R,                       # численным методом
                      z = h$Survived, X = X_mat,
                      is_aggregate = FALSE)
as_cov_est <- solve(t(lnL_R_Jac) %*% lnL_R_Jac)          # cчитаем оценку асимптотической ковариационной
                                                         # матрицы с помощью Якобианов, поскольку 
                                                         # численным методом Гессиан считается
                                                         # очень плохо
  # Реализуем тест
LM_value_1 <- t(lnL_R_grad) %*%                          # считаем статистику теста
            as_cov_est %*%                               # множителей Лагранжа
            lnL_R_grad
p_value_1 <- 1 - pchisq(LM_value_1, df = 2)              # рассчитываем p-value теста
                                                         # множителей Лагранжа

# С испольщованием регрессии на единицы

# Достанем датафрейм, содержащий
# переменные модели
d <- model.frame(model_probit)                           # все переменные

# Рассчитаем предварительные величины
z_li_est <- predict(model_probit)                                              
F_est <- pnorm(z_li_est)                              
f_est <- dnorm(z_li_est)

# Вычислим обобщенные остатки
gr <- ((d[, 1] - F_est) /                                # обобщенный остаток
       (F_est * (1 - F_est))) *
      f_est

# Считаем производные по коэффициентам
d_gamma <- apply(X_mat, 2, function(x)                   # производные по
                           {                             # регресионным коэффициентам
                             x * gr
                           })
d_t1 <- (gr * z_li_est ^ 2)                              # производная по t1
d_t2 <- (gr * z_li_est ^ 3)                              # производная по t2

# Сравним аналитические и численные производные
grad_df <- data.frame("Numeric" = lnL_R_grad,
                      "Analytical" = colSums(cbind(d_t1, 
                                             d_t2, 
                                             d_gamma)))
rownames(grad_df) <- c("t1", "t2", colnames(X_mat))
print(grad_df)

# Проводим LM тест
n <- nrow(d)                                             # число наблюдений
LM_df <- data.frame("my_ones" = rep(1, n),               # вектор из единиц 
                    "d_" = d_gamma,
                    d_t1, d_t2)                    
ones_regression <- summary(lm(my_ones~. + 0,             # регрессия на вектор единиц без константы
                           data = LM_df))       
R2 <- ones_regression$r.squared                          # коэффициент детерминации регрессии
LM_value_2 <- R2 * n                                     # LM статистика
p_value_2 <- 1 - pchisq(q = LM_value_2, 
                        df = 2)

# Сравним полученные результаты и убедимся,
# что они полностью совпадают
c(LM_value_1, LM_value_2)                                # сравниваем статистики
c(p_value_1, p_value_2)                                  # сравниваем p-value

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 1.1.    Проверьте гипотезу о нормальности при помощи:
#         1*)   LR теста
#         2*)   Wald теста
#         3**)  бутстрапа
# 1.2**.  Убедитесь, что используя Гесианный метод расчета
#         асимптотической ковариационной матрицы расчеты
#         окажутся неточными. Попытайтесь вручную
#         осуществить численное дифференцирование таким
#         образом, чтобы расчеты все же оказались точными
# 1.3***. Разработайте и реализуйте собственный тест на
#         нормальность, основанный на обобщенном
#         нормальном распределении
#         Подсказка: используйте пакет "pgnorm" и 
#                    зафиксируйте коэффициент при первом
#                    регрессоре на единице

#---------------------------------------------------
# Часть 4. Минимум Хи-квадрат метод
#---------------------------------------------------

#*******************************
# Академический уровень
#*******************************

set.seed(123)                                            # для воспроизводимости

n <- 10000                                               # число наблюдений 

# Создадим несколько переменных,
# характеризующих индивида
male <- rbinom(n, 1, 0.6)                                # переменная на пол
u_rv <- runif(n)                                         # вспомогательная переменная
higher_educ <- u_rv > 0.8                                # высшее образование
basic_educ <- (u_rv >= 0.2) & (u_rv <= 0.8)              # среднее образование
no_educ <- (u_rv < 0.2)                                  # нет образования

# Создадим сгруппированные переменные
male_higher <- male & higher_educ                        # мужчины с высшим образованием
male_basic <- male & basic_educ                          # мужчины со средним образованием
male_no <- male & no_educ                                # мужчины без образования
female_higher <- !male & higher_educ                     # женщины с высшим образованием
female_basic <- !male & basic_educ                       # женщины со средним образованием
female_no <- !male & no_educ                             # женщины без образования

# Зададим имена группам
group_names <- c("male_higher", "male_basic",
                 "male_no", "female_higher",
                 "female_basic", "female_no")

# Рассчитаем размеры групп
group_n <- c(sum(male_higher), sum(male_basic),
             sum(male_no), sum(female_higher),
             sum(female_basic), sum(female_no))

# Симулируем случайную ошибку из
# распределения стьюдента
df <- 5
eps <- rt(n, df)

# Создадим латентную переменную
gamma <- c(-2, 2, 1.5, 1)
work_li <- gamma[1] * male + gamma[2] * higher_educ +
           gamma[3] * basic_educ + gamma[4] * no_educ
work_latent <- work_li + eps
work <- as.numeric(work >= 0)

# Зададим вероятности занятости
# для каждой из групп
p_male_higher <- pt(work_li[male_higher][1], df)         # мужчины с высшим образованием
p_male_basic <- pt(work_li[male_basic][1], df)           # мужчины со средним образованием
p_male_no <- pt(work_li[male_no][1], df)                 # мужчины без образования
p_female_higher <- pt(work_li[female_higher][1], df)     # женщины с высшим образованием
p_female_basic <- pt(work_li[female_basic][1], df)       # женщины со средним образованием
p_female_no <- pt(work_li[female_no][1], df)             # женщины без образования
p <- c(p_male_higher, 
       p_male_basic,
       p_male_no,
       p_female_higher,
       p_female_basic,
       p_female_no)


# Распределим занятость в соответствии
# с заданными вероятностями
work <- rep(NA, n)
work[male_higher] <- rbinom(n = sum(male_higher),        # мужчины с высшим образованием
                            size = 1,
                            prob = p_male_higher)
work[male_basic] <- rbinom(n = sum(male_basic),          # мужчины со средним образованием
                            size = 1,
                            prob = p_male_basic)         
work[male_no] <- rbinom(n = sum(male_no),                # мужчины без образования
                           size = 1,
                           prob = p_male_no)
work[female_higher] <- rbinom(n = sum(female_higher),    # женщины с высшим образованием
                            size = 1,
                            prob = p_female_higher)
work[female_basic] <- rbinom(n = sum(female_basic),      # женщины со средним образованием 
                           size = 1,
                           prob = p_female_basic)
work[female_no] <- rbinom(n = sum(female_no),            # женщины без образования
                        size = 1,
                        prob = p_female_no)

# Оценим вероятности занятости
p_male_higher_est <- mean(work[male_higher])             # мужчины с высшим образованием
p_male_basic_est <- mean(work[male_basic])               # мужчины со средним образованием
p_male_no_est <- mean(work[male_no])                     # мужчины без образования
p_female_higher_est <- mean(work[female_higher])         # женщины с высшим образованием
p_female_basic_est <- mean(work[female_basic])           # женщины со средним образованием 
p_female_no_est <- mean(work[female_no])                 # женщины без образования

# Оценим линейное уравнение при помощи МНК
p_hat <- c(p_male_higher_est, 
           p_male_basic_est,
           p_male_no_est,
           p_female_higher_est,
           p_female_basic_est,
           p_female_no_est)
names(p_hat) <- group_names
X <- matrix(c(1, 1, 0, 0,
              1, 0, 1, 0,
              1, 0, 0, 1,
              0, 1, 0, 0,
              0, 0, 1, 0,
              0, 0, 0, 1),
            ncol = 4,
            byrow = TRUE)
colnames(X) <- c("male", "higher", "basic", "no")

# Применим обычный МНК с правильной
# квантильной функцией
model_lm <- lm(qt(p_hat, df) ~ X + 0)

# Сравним предсказанные вероятности
p_hat_lm <- pt(predict(model_lm), df)
rbind(p, p_hat, p_hat_lm)

# Сравним коэффициенты
gamma_lm <- model_lm$coefficients
data.frame("Estimates" = gamma_lm, 
           "True" = gamma)

# Применим МНК с весами

# Считаем веса
p_weight <- sqrt(p_hat * (1 - p_hat) / group_n) * 
            grad(qt, p_hat, df = df)

# Применяем МНК с весами
model_weight <- lm(qt(p_hat, df) ~ X + 0, 
                   weights = p_weight)

# Сравним предсказанные вероятности
p_hat_weight <- pt(predict(model_weight), df)
rbind(p, p_hat, p_hat_lm, p_hat_weight)

# Сравним коэффициенты
gamma_weight <- model_weight$coefficients
data.frame("LS" = gamma_lm, 
           "LS Weights" = gamma_weight,
           "True" = gamma)

# Применим двухшаговую процедуру
# Считаем веса
p_hat_weight_2 <- sqrt(p_hat_weight * 
                       (1 - p_hat_weight) / 
                       group_n) * 
                  grad(qt, p_hat_weight, df = df)

# Применяем МНК с весами
model_weight_2 <- lm(qt(p_hat, df) ~ X + 0, 
                     weights = p_hat_weight)

# Сравним предсказанные вероятности
p_hat_weight_2 <- pt(predict(model_weight_2), df)
rbind(p, p_hat, p_hat_lm, 
      p_hat_weight, p_hat_weight_2)

# Сравним коэффициенты
gamma_weight_2 <- model_weight_2$coefficients
data.frame("LS" = gamma_lm, 
           "LS Weights" = gamma_weight,
           "LS Weights 2" = gamma_weight_2,
           "True" = gamma)

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 1.1.    Рассмотрите модели с альтернативными спецификациями:
#         1)    линейная
#         2)    нормальное распределение
#         3)    логистическое распределение
#         4)    экспоненциальное распределение
# 1.2.    Примените минимум хи-квадрат метод для оценивания
#         вероятности выжить на Титанике в зависимости от пола
#         пассажира и его класса