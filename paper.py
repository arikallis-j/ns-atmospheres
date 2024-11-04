import numpy as np
from icecream import ic

A = np.array([np.float64(0.001431672094375175), np.float64(0.0030505161965432206), np.float64(0.006174836855402529), np.float64(0.0018758829493320762), np.float64(0.0018502069863216696), np.float64(0.0035562797186148055), np.float64(0.007786073867804838), np.float64(0.005489647488086043), np.float64(0.0022506402150401805), np.float64(0.002189787421277284), np.float64(0.004204739238931646), np.float64(0.002661400298658493), np.float64(0.006452709350336514), np.float64(0.006552990666310251), np.float64(0.002248795019287597), np.float64(0.005013068546543822), np.float64(0.008910979842445292), np.float64(0.0002177813650203364), np.float64(0.0017969167768819906), np.float64(0.006553846031322559), np.float64(0.008346708374590216), np.float64(0.008790687583074249), np.float64(0.00015736449614431992), np.float64(0.0016860553124143785), np.float64(0.008597188751326128), np.float64(0.0014401972098909144), np.float64(0.004087607120163954), np.float64(0.0022957210947845996), np.float64(0.009572031098679205), np.float64(0.008309122867873218), np.float64(0.00021864108905645052), np.float64(0.0012124853143275249), np.float64(0.006550934656467732), np.float64(0.0008574268768289383), np.float64(0.0018236147961109362), np.float64(0.007488012824590292), np.float64(0.0052763159818405595), np.float64(0.005791676995950222)])
X = np.arange(1.12, 50.1, 1.3)
F = np.sin(X/10.0) + A

def map1(xold, fold, xnew):
    # Определение переменных
    l = 1
    N = len(xold)

    # Определяем точку относительно массива xold
    while l < N and xnew >= xold[l]:
        l += 1
        
    # Краевой случай для l == 1 or l == N
    if l == 1 or l == N:
        ic("Border case: line interpolation")
        l = min(N - 1, l)
        x1, x0 = xold[l], xold[l-1]
        f1, f0 = fold[l], fold[l-1]
        c_l = 0.0
        b_l = (f1 - f0) / (x1 - x0)
        a_l = f1 - x1 * b_l
        return a_l + (b_l + c_l * xnew) * xnew

    # Основная квадратичная интерполяция

    x2, x1, x0 = xold[l], xold[l-1], xold[l-2]
    f2, f1, f0 = fold[l], fold[l-1], fold[l-2]

    d  = (f1 - f0) / (x1 - x0)
    c = (f2 / ((x2 - x1) * (x2 - x0))) + ((f0 / (x2 - x0)) - (f1 / (x2 - x1))) / (x1 - x0)
    b = d - (x1 + x0) * c
    a = f0 - x0 * d + x1 * x0 * c

    cm , bm, am = c, b, a

    # Дополнительная квадратичная интерполяция, если не (почти) краевой случай
    if l != N - 1:
        ic("General case: quadric interpolation with correction")
        x2, x1, x0 = xold[l+1], xold[l], xold[l-1]
        f2, f1, f0 = fold[l+1], fold[l], fold[l-1]

        d  = (f1 - f0) / (x1 - x0)
        c = (f2 / ((x2 - x1) * (x2 - x0))) + ((f0 / (x2 - x0)) - (f1 / (x2 - x1))) / (x1 - x0)
        b = d - (x1 + x0) * c
        a = f0 - x0 * d + x1 * x0 * c

        cp , bp, ap = c, b, a

        wt = abs(cp) / (abs(cp) + abs(cm)) if abs(cp) != 0 else 0.0

        a = ap + wt * (am - ap)
        b = bp + wt * (bm - bp)
        c = cp + wt * (cm - cp)

    else:
        ic("Specific case: quadric interpolation without correction")

    return a + (b + c * xnew) * xnew

print(F)
print(X)
R = map1(X, F, 53.4)
print(R)
R = map1(X, F, 1.14)
print(R)
R = map1(X, F, 49.21)
print(R)
R = map1(X, F, 2.43)
print(R)




# def map1(x_old, f_old, n_old, x_new):
#     # Инициализируем переменные и размеры
#     l = 2
#     ll = 0
#     a, b, c = 0.0, 0.0, 0.0
    
#     # Найдем подходящий индекс l
#     while l <= n_old and x_new >= x_old[l - 1]:
#         l += 1

#     # Обработка случая выхода за пределы массива
#     if l > n_old:
#         if l == ll:
#             return a + (b + c * x_new) * x_new, ll - 1
#         l = min(n_old, l)
#         c = 0.0
#         b = (f_old[l - 1] - f_old[l - 2]) / (x_old[l - 1] - x_old[l - 2])
#         a = f_old[l - 1] - x_old[l - 1] * b
#         ll = l
#         return a + (b + c * x_new) * x_new, ll - 1

#     # Если ll совпадает с l, используем текущие значения
#     if l == ll:
#         return a + (b + c * x_new) * x_new, ll - 1

#     # Проверка для l == 2 (особый случай)
#     if l == 2:
#         l = min(n_old, l)
#         c = 0.0
#         b = (f_old[l - 1] - f_old[l - 2]) / (x_old[l - 1] - x_old[l - 2])
#         a = f_old[l - 1] - x_old[l - 1] * b
#         ll = l
#         return a + (b + c * x_new) * x_new, ll - 1

#     # Основная интерполяция
#     l1 = l - 1
#     if l > ll + 1 or l == 3:
#         l2 = l - 2
#         d = (f_old[l1 - 1] - f_old[l2 - 1]) / (x_old[l1 - 1] - x_old[l2 - 1])
#         cbac = f_old[l - 1] / ((x_old[l - 1] - x_old[l1 - 1]) * (x_old[l - 1] - x_old[l2 - 1])) + \
#                (f_old[l2 - 1] / (x_old[l - 1] - x_old[l2 - 1]) - f_old[l1 - 1] / (x_old[l - 1] - x_old[l1 - 1])) / \
#                (x_old[l1 - 1] - x_old[l2 - 1])
#         bbac = d - (x_old[l1 - 1] + x_old[l2 - 1]) * cbac
#         abac = f_old[l2 - 1] - x_old[l2 - 1] * d + x_old[l1 - 1] * x_old[l2 - 1] * cbac
#         if l == n_old:
#             c = cbac
#             b = bbac
#             a = abac
#             ll = l
#             return a + (b + c * x_new) * x_new, ll - 1

#     # Если l < n_old, используем расчет на основе cfor, bfor, afor
#     d = (f_old[l - 1] - f_old[l1 - 1]) / (x_old[l - 1] - x_old[l1 - 1])
#     cfor = f_old[l] / ((x_old[l] - x_old[l - 1]) * (x_old[l] - x_old[l1 - 1])) + \
#            (f_old[l1 - 1] / (x_old[l] - x_old[l1 - 1]) - f_old[l - 1] / (x_old[l] - x_old[l - 1])) / \
#            (x_old[l - 1] - x_old[l1 - 1])
#     bfor = d - (x_old[l - 1] + x_old[l1 - 1]) * cfor
#     afor = f_old[l1 - 1] - x_old[l1 - 1] * d + x_old[l - 1] * x_old[l1 - 1] * cfor
#     wt = 0.0
#     if abs(cfor) != 0.0:
#         wt = abs(cfor) / (abs(cfor) + abs(cbac))
#     a = afor + wt * (abac - afor)
#     b = bfor + wt * (bbac - bfor)
#     c = cfor + wt * (cbac - cfor)
#     ll = l
#     return a + (b + c * x_new) * x_new, ll - 1


def map2(X, F, Y):
    # Приводим массивы к нужным формам для векторных вычислений
    X = np.array(X)  # X имеет размер (L,)
    F = np.array(F)  # F имеет размер (L, P, A, B)
    Y = np.array(Y)  # Y имеет размер (A, B)
    
    # Получаем размеры массивов
    l_shape = X.shape[0]
    p_shape = F.shape[1]
    g_shape_a, g_shape_b = Y.shape

    # Векторизуем X и F для соответствия размерности (L, P * A * B)
    X = np.tile(X[:, None], (1, p_shape * g_shape_a * g_shape_b))
    F = F.reshape((l_shape, p_shape * g_shape_a * g_shape_b))
    Y = Y.flatten()  # Преобразуем Y в одномерный массив для удобства
    
    # Инициализация для результата
    R = np.zeros_like(Y, dtype=float)  # Итоговый результат, который будем заполнять
    
    # Вспомогательные массивы для коэффициентов
    a = np.zeros_like(Y, dtype=float)
    b = np.zeros_like(Y, dtype=float)
    c = np.zeros_like(Y, dtype=float)

    # Создаем массивы для хранения промежуточных данных
    L = np.searchsorted(X[:, 0], Y, side='right')  # Индексы интерполяции для каждого Y
    
    # Логические маски для разных случаев
    mask1 = (L >= l_shape) | (L < 2)
    mask2 = (L == 2) | (L == l_shape - 1)
    mask3 = ~mask1 & ~mask2  # Общий случай для квадратичной интерполяции

    # Граничный случай: линейная интерполяция для крайних точек
    a[mask1] = F[L[mask1] - 1]
    b[mask1] = (F[L[mask1]] - F[L[mask1] - 1]) / (X[L[mask1]] - X[L[mask1] - 1])
    c[mask1] = 0.0


    # Квадратичная интерполяция для основного случая
    d = (F[L - 1] - F[L - 2]) / (X[L - 1] - X[L - 2])
    cbac = F[L] / ((X[L] - X[L - 1]) * (X[L] - X[L - 2])) + \
           (F[L - 2] / (X[L] - X[L - 2]) - F[L - 1] / (X[L] - X[L - 1])) / (X[L - 1] - X[L - 2])
    bbac = d - (X[L - 1] + X[L - 2]) * cbac
    abac = F[L - 2] - X[L - 2] * d + X[L - 1] * X[L - 2] * cbac

    # Вторая часть: интерполяция c использованием ближайших значений
    afor = F[L - 1] - X[L - 1] * d + X[L] * X[L - 1] * cbac
    bfor = d - (X[L] + X[L - 1]) * cbac
    cfor = F[L + 1] / ((X[L + 1] - X[L]) * (X[L + 1] - X[L - 1])) + \
           (F[L - 1] / (X[L + 1] - X[L - 1]) - F[L] / (X[L + 1] - X[L])) / (X[L] - X[L - 1])

    # Расчет весов и итоговых коэффициентов
    wt = np.abs(cfor) / (np.abs(cfor) + np.abs(cbac))
    a[mask3] = afor + wt * (abac - afor)
    b[mask3] = bfor + wt * (bbac - bfor)
    c[mask3] = cfor + wt * (cbac - cfor)

    # Итоговая интерполяция
    R = a + b * Y + c * Y**2
    
    # Возвращаем результат в нужной форме
    return R.reshape((p_shape, g_shape_a, g_shape_b))

