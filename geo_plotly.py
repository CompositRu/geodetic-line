import plotly.graph_objs as go
import numpy as np
from math import sin, cos, asin, radians, degrees, tan, sqrt
import math


# -------------------------------------------------------------------------------------------
# ------- НАХОЖДЕНИЕ КООРДИНАТ ГЕОДЕЗИЧЕСКОЙ ЛИНИИ НА ПОВЕРХНОСТИ ТЕЛА ВРАЩЕНИЯ -------------
# -------------------------------------------------------------------------------------------


def get_shift_angle(np_fi, num_turns):
    central_angle = 2 * np.pi / num_turns
    delta_fi = (np_fi[-1] - np_fi[0]) % (2 * np.pi) 
    num_central_angles = int(delta_fi // central_angle) + 1
	# Для равномерного покрытия 
    while (math.gcd(num_turns, num_central_angles) != 1 or num_turns == num_central_angles):
        num_central_angles += 1
    shift_angle = central_angle * num_central_angles - delta_fi
    print('shift_angle rad:', shift_angle)
    print('shift_angle deg:', degrees(shift_angle))
    return shift_angle


def get_rare(numpy_ar):

    def f(ar):
        idx1 = len(ar) // 1000
        idx2 = len(ar) // 100
        idx3 = len(ar) // 10
        ar1 = ar[:idx1]
        ar2 = ar[idx1:idx2:3]
        ar3 = ar[idx2:idx3:30]
        ar4 = ar[idx3::300]
        if type(ar) == 'numpy.ndarray':
            return np.concatenate((ar1, ar2, ar3, ar4))
        else:
            return ar1 + ar2 + ar3 + ar4

    def rev(ar):
        if type(ar) == 'numpy.ndarray':
            return np.flip(ar)
        else:
            return list(reversed(ar))
  
    idx_ar = list(range(len(numpy_ar )))        
    # print(f(idx_ar[:len(idx_ar)//2]))
    # print(rev(f(rev(idx_ar[len(idx_ar)//2:]))))
    idx = f(idx_ar[:len(idx_ar)//2]) + rev(f(rev(idx_ar[len(idx_ar)//2:])))
    # print('#', len(numpy_ar[:idx1]))
    # print('#', len(numpy_ar[idx1:idx2:3]))
    # print('#', len(numpy_ar[idx2:idx3:30]))
    # print('#', len(numpy_ar[idx3::300]))
    return numpy_ar[idx]

# Начальные геометрические условия
start_z_data = [15, 80, 118, 258, 308]
start_z_data = [x - 15 for x in start_z_data]
start_diameters_data = [91, 150, 170, 236, 236]
start_r_data = [diam/2 for diam in start_diameters_data]
#z_data = [0, 15, 80, 118, 258, 308]
#diameters_data = [91, 91, 150, 170, 236, 236]
#z_data = [15, 80, 118, 258, 286, 290, 292, 294, 308]
#diameters_data = [91, 150, 170, 236, 106, 100, 100, 106, 236]

# Зеркально отображаем данные и дополняем стартовые массивы
z_reverse = [start_z_data[-1] * 2 - x for x in reversed(start_z_data[:-1])]
start_z_data = start_z_data + z_reverse
start_r_data = start_r_data + list(reversed(start_r_data[:-1]))

# Линейная интерполяция для получения num_points точек
num_points = (start_z_data[-1] - start_z_data[0] + 1) * 10
np_z = np.linspace(start_z_data[0], start_z_data[-1], num_points)
np_r = np.interp(np_z, start_z_data, start_r_data)

# Дифференциируем массив радиусов
diff_r_points = np.diff(np_r)
diff_r_points = np.append(diff_r_points, diff_r_points[-1])

# Находим коэффициент альфа, зависящий от угла армирования
# Угол армирования psi_0 задаём для конкретного радиуса r_0
psi_0 =89
r_0 = np_r[0]
alfa = 1 / (r_0 * math.sin(psi_0 * np.pi / 180) )

# Находим dz для интеграла
dz = (start_z_data[-1] - start_z_data[0]) / num_points

# Вычисляем подыинтегральное выражение (65) со стр. 426 книги Бухгольца
# для каждой точки
f = np_r
dfdz = diff_r_points
geo = dz / f * np.sqrt((1 + dfdz * dfdz) / (alfa * alfa * f * f - 1))

# Находим углы как интеграл в каждой точке
np_fi = np.cumsum(geo)

# Опрелеляем угол армировния в большем сечении через формулу Клеро 
# r * sin = const
r_center = max(np_r)
psi_center = asin(r_0 / r_center * sin(radians(psi_0)))

# Выбор количества витков    
num_turns = 50 #количество витков для полного покрытия наибольшего диаметра (взято с потолка)
width = 2 * np.pi * r_center / num_turns 
print('width:', width) 
width_center = width / cos(psi_center) # данное значение должно быть близко к реальному значению ширины жгута
print('width_center:', width_center)

# Находим смещение второго витка, чтобы все витки равномерно покрывали поверхность
shift_angle = get_shift_angle(np_fi, num_turns)

'''
# Добавляем второй виток
# Добавляем витки
np_r = np_r_1_turn
np_fi = np_fi_1_turn
np_z = np_z_1_turn

np_fi_2 = np_fi_1_turn + np_fi[-1] + shift_angle
np_r_2_turns = np.concatenate((np_r_1_turn, np_r_1_turn))
np_fi_2_turns = np.concatenate((np_fi_1_turn, np_fi_2))
np_z_2_turns = np.concatenate((np_z_1_turn, np.flip(np_z_1_turn)))

# Перед разряжение массивов данных создаём отдельно массив данных для будущего рассчёта векторов схода нити
# Объединяем r, fi, z в один массив данных и добавляем в каждую точку координаты следующей (чтобы эти данные не потерялись при расзряжении)
s = np.stack((np_r, np_fi, np_z), axis=1)
print(s)
r = np.stack((np.roll(np_r, -1), np.roll(np_fi, -1), np.roll(np_z, -1)), axis=1)
print(r)
c = np.concatenate((s, r), axis = 1) 

directions = 
'''

# Делаем массив точек разреженным
print('old len(np_fi)', len(np_fi))
np_r_1_turn = get_rare(np_r)
np_fi_1_turn = get_rare(np_fi)
np_z_1_turn = get_rare(np_z)
print('new len(np_fi)', len(np_fi_1_turn))

# Добавляем витки
np_r = np_r_1_turn
np_fi = np_fi_1_turn
np_z = np_z_1_turn

for i in range(int(num_turns - 1)):
    np_fi_2 = np_fi_1_turn + np_fi[-1] + shift_angle
    np_r = np.concatenate((np_r, np_r_1_turn))
    np_fi = np.concatenate((np_fi, np_fi_2))
    if i % 2 == 0:
        np_z = np.concatenate((np_z, np.flip(np_z_1_turn)))
    else:
        np_z = np.concatenate((np_z, np_z_1_turn))

print('len(np_fi)', len(np_fi))
print('len(np_r)', len(np_r))
print('len(np_z)', len(np_z))


# -------------------------------------------------------------------------------------------
# ------- ПОСТРОЕНИЕ 3D МОДЕЛИ --------------------------------------------------------------
# -------------------------------------------------------------------------------------------


def rotation_body(radius_array, z_array, nt = 50):
    theta = np.linspace(0, 2*np.pi, nt)
    np_fi, r = np.meshgrid(theta, radius_array)
    _,  z = np.meshgrid(theta, z_array)
    x = r * np.cos(np_fi)
    y = r * np.sin(np_fi)
    z = z
    return x, y, z    


def cartesian_from_polar(r, np_fi, h):
    x = r * np.cos(np_fi)
    y = r * np.sin(np_fi)
    z = h
    return x, y, z

def cartesian_from_polar_2(r, phi):
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return(x, y)


# Точки для построения поверхности вращения
x1, y1, z1 = rotation_body(np_r_1_turn, np_z_1_turn)
surface = go.Surface(x=z1, y=y1, z=x1, colorscale='Viridis', opacity=0.5)

# Точки для построения геодезической линии
x2, y2, z2 = cartesian_from_polar(np_r, np_fi, np_z)
curve = go.Scatter3d(x=z2, y=y2, z=x2, mode='lines', line={'width': 10})

# point = np.array([118., 3.11923172, 243.02675072])
# direct = np.array([0., 0.0003535, 0.09984665])
# x3, y3, z3 = cartesian_from_polar(point[0], point[1], point[2])
# P = go.Scatter3d(x=np.array([z3]), y=np.array([y3]), z=np.array([x3]))
# x4, y4, z4 = cartesian_from_polar(direct[0], point[1] + direct[1], direct[2])
# k = 20;
# D = go.Scatter3d(x=np.array([z4 * k + z3]) , y=np.array([y4 * k + y3]), z=np.array([x4 * k + x3]))

# Рисуем графики
# data = [surface, curve, P, D]
data = [surface, curve]
fig = go.Figure(data=data)
fig.update_layout(
    title={
        'text': f'Углы армирования:<br>на краю - {psi_0}<br>в центре - {round(degrees(psi_center))}',
        'y':0.9,
        'x':0.5,
        'xanchor': 'center',
        'yanchor': 'top'})
fig.write_html('tmp.html', auto_open=True) # Для получения графиков убираем # с этой строки


'''
# -------------------------------------------------------------------------------------------
# ------- ПОЛУЧАЕМ МАССИВ КООРДИНАТ УПРАВЛЯЮЩИХ ОРГАНОВ СТАНКА ------------------------------
# -------------------------------------------------------------------------------------------

z_plane = -175 # высота плоскости, в которой движется кольцо раздаточной головы


def cartesian_from_polar_for_one_point(r, fi, h):
    x = r * math.cos(fi)
    y = r * math.sin(fi)
    z = h
    return x.tolist(), y.tolist(), z


def intersect_3D(point, ray, point_on_plain, normal):
    t = np.dot((point_on_plain - point), normal) / np.dot(ray, normal)
    return point + t * ray


def intersect_3D_horizontal_plane(point, ray, z_plain):
    t = (z_plain - point[2]) / ray[2]
    return point + t * ray


def get_distance_on_horizontal_plane(a, b, point):
    ax, ay, bx, by, x, y = *a, *b, *point
    a = ay - by
    b = bx - ax
    c = ax * by - bx * ay
    return abs(a * x + b * y + c) / sqrt(a * a + b * b)


def get_distace(segments, point):
    distances  = []
    for i in range(len(segments) - 1):
        line = (segments[i], segments[i + 1])
        distances.append(get_distance_on_horizontal_plane(*line, point))
    return min(distances)


def get_vector_direction(i):
    i = i % ( 2 * len(np_r_1_turn))
    direction = 'forward' if i < len(np_r_1_turn) else 'backward'

    r0 = np_r_1_turn[i]
    fi0 = np_fi_1_turn[i]
    r1 = np_r_1_turn[i + 1] if i != len(np_r_1_turn) - 1 else np_r_1_turn[0]
    fi1 = np_fi_1_turn[i + 1] if i != len(np_r_1_turn) - 1 else np_fi_1_turn[0]

    i = i if direction == 'forward' else len(np_r_1_turn) - i
    next_i = i + 1 if direction == 'forward' else i - 1
    z0 = np_z_1_turn[i] 
    z1 = np_z_1_turn[next_i] if i != len(np_r_1_turn) - 1 else np_z_1_turn[0]      

# Получаем массив точек без разряжения   
np_r = np_r_1_turn
np_fi = np_fi_1_turn
np_z = np_z_1_turn
print(len(np_r))
print(len(np_fi))
print(len(np_z))

for i in range(int(num_turns - 1)):
    np_fi_2 = np_fi_1_turn + np_fi[-1] + shift_angle
    np_r = np.concatenate((np_r, r))
    np_fi = np.concatenate((np_fi, np_fi_2))
    if i % 2 == 0:
        np_z = np.concatenate((np_z, np.flip(z)))
    else:
        np_z = np.concatenate((np_z, z))

print('len(np_r)', len(np_r))
print('len(np_fi)', len(np_fi))
print('len(np_z)', len(np_z))

# x, y, z = cartesian_from_polar(np_r, np_fi, np_z)
# points = np.stack((x, y, z)).T
points = np.stack((np_r, np_fi, np_z)).T # либо вместо траспонирования ставим параметр axis=-1

# Находим вектора направления схода нити в каждой точке относительно предыдущей
directions = np.diff(points, axis=0)
directions = np.append(directions, [directions[-1]], axis = 0)

# Траектрия движения раздаточной головы. Условно считаем, что идеальная траектория - эквидестанта поверхности.
# Эквидестанту получаем простым параллельным переносом
# По краям эквидестанты добавляем две точки для увеличения длины траектории
shift = 50 # сдвиг
increasing_distance = 100 # Увеличиваем вылет траектории по краям
first_point = [(-increasing_distance, start_r_data[0] + shift)] 
curve = [(x, y + shift) for x, y in zip(start_z_data, start_r_data)]
end_point = [(start_z_data[-1] + increasing_distance, start_r_data[-1] + shift)]
trajectory = first_point + curve + end_point
print(trajectory)

max_distance = 0; # параметр исключительно для проверки результатов
# for i in range(len(points)):
for i in range(2 * len(np_r_1_turn)):
    r = points[i][0]
    fi = points[i][1]
    z = points[i][2]
    dr = directions[i][0]
    dfi = directions[i][1]
    dz = directions[i][2]

    start_beta = radians(-80)
    finish_beta = radians(45)
    steps = 200
    dbeta = (finish_beta - start_beta) / steps;
    min_distance = None
    if i == 2434:
        print(points[i], directions[i])
        break
    # for step in range(steps):
    #     beta = start_beta + step * dbeta;
    #     y, z, x = cartesian_from_polar(r, beta, z) # оси координат совпадают с осями зелёного станка
    #     ry, rz, rx = cartesian_from_polar(dr, beta + dfi, dz)
    #     if abs(rz) < float(0.0001): 
    #         # print(i, step, rz)
    #         continue
    #     p = intersect_3D_horizontal_plane(np.array([x, y, z]), np.array([rx, ry, rz]), z_plane).tolist()
    #     d = get_distace(trajectory, (p[0], p[1]))
    #     if step == 0:
    #         min_distance = d
    #     else:
    #         min_distance = d if d < min_distance else min_distance
    # max_distance = min_distance if min_distance > max_distance else max_distance
    # if i % 1000 == 0: print(i, max_distance)
print(max_distance)
'''