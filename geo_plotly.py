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
    return shift_angle


def get_rare(ar):
    # Разряжаем массив опираясь на внутреннее чутьё
    def f(ar):
        idx1 = len(ar) // 1000
        idx2 = len(ar) // 100
        idx3 = len(ar) // 10
        ar1 = ar[:idx1:5]
        ar2 = ar[idx1:idx2:30]
        ar3 = ar[idx2:idx3:100]
        ar4 = ar[idx3::350]
        if type(ar) == 'numpy.ndarray':
            return np.concatenate((ar1, ar2, ar3, ar4))
        else:
            return ar1 + ar2 + ar3 + ar4
    def rev(ar):
        if type(ar) == 'numpy.ndarray':
            return np.flip(ar)
        else:
            return list(reversed(ar)) 
    idx_ar = list(range(len(ar)))        
    idx = f(idx_ar[:len(idx_ar)//2]) + rev(f(rev(idx_ar[len(idx_ar)//2:])))
    if type(ar) == 'numpy.ndarray':
        return numpy_ar[idx]
    else:
        res = []
        for i in idx:
            res.append(ar[i])
        return res


# Начальные геометрические условия

start_z_data = [15, 80, 118, 258, 308] # Стартовые данные
start_z_data = [x - 15 for x in start_z_data] # Стартовые данные
start_diameters_data = [91, 150, 170, 236, 236] # Стартовые данные
start_r_data = [diam / 2 for diam in start_diameters_data] # Стартовые данные

# Зеркально отображаем данные и дополняем стартовые массивы
z_reverse = [start_z_data[-1] * 2 - x for x in reversed(start_z_data[:-1])]
start_z_data = start_z_data + z_reverse
start_r_data = start_r_data + list(reversed(start_r_data[:-1]))

# Задаём угол армирования на краях. В идеале он должен быть равен 90, но в этом случае слишком много нитки наматывается на краях
psi_0 = 89

def get_geodetic(r_data, z_data):
    # Линейная интерполяция для получения num_points точек
    num_points = (z_data[-1] - z_data[0] + 1) * 10
    np_z = np.linspace(z_data[0], z_data[-1], num_points)
    np_r = np.interp(np_z, z_data, r_data)

    # Дифференциируем массив радиусов
    diff_r_points = np.diff(np_r)
    diff_r_points = np.append(diff_r_points, diff_r_points[-1])

    # Находим коэффициент альфа, зависящий от угла армирования
    # Угол армирования psi_0 задаём для конкретного радиуса r_0
    r_0 = np_r[0]
    alfa = 1 / (r_0 * math.sin(radians(psi_0)))

    # Находим dz для интеграла
    dz = (z_data[-1] - z_data[0]) / num_points

    # Вычисляем подыинтегральное выражение (65) со стр. 426 книги Бухгольца
    # для каждой точки
    f = np_r
    dfdz = diff_r_points
    geo = dz / f * np.sqrt((1 + dfdz * dfdz) / (alfa * alfa * f * f - 1))

    # Находим углы как интеграл в каждой точке
    np_fi = np.cumsum(geo)
    np_fi -= np_fi[0]

    return np_r, np_fi, np_z

np_r, np_fi, np_z = get_geodetic(start_r_data, start_z_data)

# Опрелеляем угол армировния в большем сечении через формулу Клеро 
# r * sin = const
r_center = max(np_r)
r_0 = np_r[0]
psi_center = asin(r_0 / r_center * sin(radians(psi_0)))

# Выбор количества витков    
num_turns = 100 #количество витков для полного покрытия наибольшего диаметра (взято с потолка)
width = 2 * np.pi * r_center / num_turns 
print('width:', width) 
width_center = width / cos(psi_center) # данное значение должно быть близко к реальному значению ширины жгута
print('width_center:', width_center)

# Находим смещение второго витка, чтобы все витки равномерно покрывали поверхность
shift_angle = get_shift_angle(np_fi, num_turns)

# Добавляем второй виток
np_r_1_turn = np_r
np_fi_1_turn = np_fi
np_z_1_turn = np_z

np_fi_2 = np_fi_1_turn + np_fi[-1] + shift_angle
np_r_2_turns = np.concatenate((np_r_1_turn, np_r_1_turn))
np_fi_2_turns = np.concatenate((np_fi_1_turn, np_fi_2))
np_z_2_turns = np.concatenate((np_z_1_turn, np.flip(np_z_1_turn)))

# Перед разряжение массивов данных создаём отдельно массив данных для будущего рассчёта векторов схода нити
# Объединяем r, fi, z в один массив данных и добавляем в каждую точку координаты следующей (чтобы эти данные не потерялись при расзряжении)
directions_1 = np.stack((np_r_2_turns, np_fi_2_turns, np_z_2_turns), axis=1) # либо axis=1, либо транспонируем полученную матрицу 
directions_2 = np.stack((np.roll(np_r_2_turns, -1), np.roll(np_fi_2_turns, -1), np.roll(np_z_2_turns, -1)), axis=1)
directions = np.concatenate((directions_1, directions_2), axis = 1)

# Т.к. значение полярного угла для точек всегда растёт (т.к. вычисляется как интеграл),
# то неверно использовать для направления касательной в последней точке координату полярного угла в первой точке.
# Хотя практика показала, что даже большая погрешность значения в этой ячейке массива практически не влияет на результат,
# проводим элементарную замену значения в массиве.
directions[len(directions) - 1][4] = directions[len(directions) - 1][1] + shift_angle

# Для разряжения массива np_array подготовим список индексов
range_1 = list(range(len(np_r_1_turn))) 
range_2 = list(range(len(np_r_1_turn), 2 * len(np_r_1_turn))) 
range_1 = get_rare(range_1)
range_2 = get_rare(range_2)
idx_range = range_1 + range_2

# Делаем массив точек разреженным
print('old len(np_r_2_turn)', len(np_r_2_turns))
np_r_2_turns = np_r_2_turns[idx_range]
np_fi_2_turns = np_fi_2_turns[idx_range]
np_z_2_turns = np_z_2_turns[idx_range]
directions = directions[idx_range]
print('new len(np_r_2_turn)', len(np_r_2_turns))

# Создаём массивы из num_turns витков
np_r = np.array([])
np_fi = np.array([])
np_z = np.array([])

double_turns = int(num_turns / 2)
for i in range(double_turns):    
    np_fi_2 = np_fi_2_turns if i == 0 else np_fi_2_turns + np_fi[-1] + shift_angle
    np_r = np.concatenate((np_r, np_r_2_turns))
    np_fi = np.concatenate((np_fi, np_fi_2))
    np_z = np.concatenate((np_z, np_z_2_turns))

# -------------------------------------------------------------------------------------------
# ------- ПОСТРОЕНИЕ ТРАЕКТОРИИ ДВИЖЕНИЯ РАСКЛАДОЧНОЙ ГОЛОВЫ --------------------------------
# -------------------------------------------------------------------------------------------

# Траектрия движения раздаточной головы. Условно считаем, что идеальная траектория - эквидестанта поверхности.
# Эквидестанту получаем простым параллельным переносом
# По краям эквидестанты добавляем две точки для увеличения длины траектории
z_plane = -140 # высота плоскости, в которой движется кольцо раздаточной головы
shift = 0 # сдвиг
increasing_distance = 154 # Увеличиваем вылет траектории по краям
extreme_y = 0
first_point = [(-increasing_distance, extreme_y)] 
curve_trajectory = [[x, y + shift] for x, y in zip(start_z_data, start_r_data)]
end_point = [(start_z_data[-1] + increasing_distance, extreme_y)]
trajectory_list = first_point + curve_trajectory + end_point

# Итоговая траектория близка к эквидестанте и подобрана вручную, опираясь на внутреннее чутьё
trajectory_list[1][1] += 25
trajectory_list[-2][1] += 25
trajectory_list[2][1] += 15
trajectory_list[-3][1] += 15
trajectory_list[3][1] += 10
trajectory_list[-4][1] += 10

# -------------------------------------------------------------------------------------------
# ------- ПОСТРОЕНИЕ 3D МОДЕЛИ И ПОЛУЧЕНИЕ КООРДИНАТ УПРАВЛЯЮЩИХ ОРГАНОВ СТАНКА -------------
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


def inside(x, interval):
    if x >= interval[0] and x <= interval[1]:
        return True
    return False


def intersection_plane(point, ray, z_plain):
    t = (z_plain - point[2]) / ray[2]
    if t < 0:
        print('t < 0', point) # Необходимо учитывать, что при t< 0 мы получаем некореектные данные (вектор в противоположную сторону направлен) и их надо отбрасывать
    return point + t * ray


def get_distance_on_horizontal_plane(a, b, point):
    np_a = np.array(a)
    np_b = np.array(b)
    np_p = np.array(point) if type(point) != 'numpy.ndarray' else point
    np_p = np_p[:2] if len(np_p) > 2 else np_p
    length = np.linalg.norm
    if np.dot(np_b - np_a, np_p - np_a) > 0:
        if np.dot(np_a - np_b, np_p - np_b) > 0:
            dx, dy = np_a - np_b
            c = np_a[0] * np_b[1] - np_b[0] * np_a[1]
            return abs(dy * np_p[0] - dx * np_p[1] + c) / length(np_b - np_a)
        else:
            return length(np_p - np_b)
    else:
        return length(np_p - np_a)


def get_distance(segments, point):
    distances  = []
    for i in range(len(segments) - 1):
        line = (segments[i], segments[i + 1])
        distances.append(get_distance_on_horizontal_plane(*line, point))
    return min(distances)


def get_intersect_point_exit_point_angle(start_polar, end_polar):
    dfi = (end_polar - start_polar)[1]
    start_beta = radians(45)
    finish_beta = radians(150)
    steps = 200
    dbeta = (finish_beta - start_beta) / steps;
    min_distance = None
    point_on_plain = None
    start_point = None
    angle = start_beta
    for step in range(steps):
        beta = start_beta + step * dbeta;
        z1, y1, x1 = cartesian_from_polar(start_polar[0], beta, start_polar[2])
        z2, y2, x2 = cartesian_from_polar(end_polar[0], beta + dfi, end_polar[2])
        np_start = np.array((x1, y1, z1))
        np_end = np.array((x2, y2, z2))
        np_ray = np_end - np_start
        intersection_point = intersection_plane(np_start, np_ray, z_plane)
        d = get_distance(trajectory_list, intersection_point)
        if step == 0:
            min_distance = d
            point_on_plain = intersection_point
            start_point = np_start
        else:
            if d < min_distance:
                min_distance = d
                point_on_plain = intersection_point
                start_point = np_start
                angle = beta
    return point_on_plain, start_point, angle, min_distance


# Точки для построения поверхности вращения
x1, y1, z1 = rotation_body(np_r_1_turn, np_z_1_turn)
surface = go.Surface(x=z1, y=y1, z=x1, colorscale='Viridis', opacity=0.5, name="surface")

# Точки для построения геодезической линии
x2, y2, z2 = cartesian_from_polar(np_r, np_fi, np_z)
curve = go.Scatter3d(visible='legendonly', x=z2, y=y2, z=x2, mode='lines', line={'width': 10}, name="curve")

# Точки для построения траектории движения раскладочной головы
np_trajectory = np.array(trajectory_list)
z_column = np.array([z_plane] * len(np_trajectory))
np_trajectory = np.column_stack((np_trajectory, z_column)).T
x, y, z = np_trajectory[0], np_trajectory[1], np_trajectory[2]
trajectory = go.Scatter3d(x=x, y=y, z=z, mode='lines', line={'width': 10}, name="trajectory")

# Точки для построения границ трактории
center = (np.max(np_trajectory[0]) - np.min(np_trajectory[0])) / 2 + np.min(np_trajectory[0])
work_area = 900 # Длина рабочей зоны раскладочной головы до срабатывания концевиков
x1 = center - work_area / 2
x = np.array([x1, x1])
y = np.array([-np.max(np_r_1_turn), np.max(np_r_1_turn)])
z = np.array([z_plane, z_plane])
border_1 = go.Scatter3d(x=x, y=y, z=z, mode='lines', line={'width': 10}, marker={'color':'red'}, name="borders", legendgroup="borders", showlegend=True)
x2 = center + work_area / 2
x = np.array([x2, x2])
border_2 = go.Scatter3d(x=x, y=y, z=z, mode='lines', line={'width': 10}, marker={'color':'red'}, name="borders", legendgroup="borders", showlegend=False)
borders = [border_1, border_2]
borders_interval = [x1.tolist(), x2.tolist()]
print('borders', [x1, x2])

# Точки для построения векторов схода нити
vectors = [] # Массив для построения линий схода нити в 3D
thread_points = [] # Массив некоторых координат линий схода нити для последующего получения координат управлющих органов станка
exit_line = [] # Массив точек схода нити
max_distance = 0; # Параметр исключительно для проверки результатов
for i, item in enumerate(directions):
    start_polar = item[:3]
    end_polar = item[3:]
    point_on_plain, start_point, angle, min_distance = get_intersect_point_exit_point_angle(start_polar, end_polar)
    x, y, fi = point_on_plain[0].tolist(), point_on_plain[1].tolist(), angle
    thread_points.append([x, y, fi])
    exit_line.append(start_point)
    v = np.stack((start_point, point_on_plain)).T
    x, y, z = v[0], v[1], v[2]
    color = 'cyan' if inside(point_on_plain[0].tolist(), borders_interval) else 'coral'
    flag = True if i == 0 else False
    vectors.append(go.Scatter3d(x=x, y=y, z=z, mode='lines', line={'width': 5}, opacity=0.3, marker={'color':color}, legendgroup="group", showlegend=flag, name="threads")) # showlegend=False
    max_distance = min_distance if min_distance > max_distance else max_distance
print('max_distance', max_distance) # параметр исключительно для проверки

# Точки для построения для линии схода нити
x, y, z = np.array(exit_line).T
exit_line = go.Scatter3d(x=x, y=y, z=z, mode='lines', line={'width': 10}, opacity=0.8, marker={'color':'#FF8C00'}, name="exit_line")

# Создание html файла 3D модели
data = [surface, curve, trajectory, *borders, *vectors, exit_line]
fig = go.Figure(data=data)
fig.update_layout(legend_orientation="h", 
                  legend=dict(x=.5, xanchor="center"),
                  margin=dict(l=0, r=0, t=0, b=0))
fig.update_layout(
    title={ 
        'text': f'Углы армирования:<br>на краю - {psi_0}<br>в центре - {round(degrees(psi_center))}',
        'y':0.9,
        'x':0.5,
        'xanchor': 'center',
        'yanchor': 'top'})
fig.write_html('tmp.html', auto_open=True) # Для получения графиков убираем # с этой строки


# -------------------------------------------------------------------------------------------
# ------- ЗАПИСЬ В ФАЙЛ МАССИВА КООРДИНАТ УПРАВЛЯЮЩИХ ОРГАНОВ СТАНКА ------------------------
# -------------------------------------------------------------------------------------------

shift_Y_axis = 100 # Смещение реальной СК, относительно той, в которой проводился расчёт
shift_X_axis = (np.max(np_z_1_turn) - np.min(np_z_1_turn)) / 2 + np.min(np_z_1_turn) # На станке выставляем ноль по центру оснастки
print('length', shift_X_axis * 2)
print('np.min(np_z_1_turn)', np.min(np_z_1_turn))
print('np.max(np_z_1_turn)', np.max(np_z_1_turn))
borders_interval = [x - shift_X_axis for x in borders_interval]

with open("Winding_frame.JOB", "w") as f:
     # Аккуратный подъезд на стартовое положение
    n = 804
    f.write(f'n{n} G01 X{round(thread_points[0][0] - shift_X_axis , 1)} Y{round(thread_points[0][1] + shift_Y_axis + 150, 1)}\n')
    n += 1
    f.write(f'n{n} G01 Y{round(thread_points[0][1] + shift_Y_axis, 1)}\n')
    n += 1
    beta_start = 30.0 + 20158.8
    beta = beta_start
     # Поворот на beta_start для укладки нити в своё стартовое положение
    f.write(f'n{n} G01 A{-beta}\n')
    n += 1
    real_fi_prev = None
    curve_fi_prev = None

    for turn in range(double_turns):
        for i, item in enumerate(thread_points):
            if turn == 0 and i == 0:
                real_fi_prev = degrees(item[2])
                curve_fi_prev = degrees(np_fi[0])
            else:
                real_fi_prev = real_fi
                curve_fi_prev = curve_fi
            real_fi = degrees(item[2])
            curve_fi = degrees(np_fi[i + turn * len(np_r_2_turns)])
            dbeta = (curve_fi - curve_fi_prev) - (real_fi - real_fi_prev)
            beta += dbeta
            x = round(item[0] - shift_X_axis , 1)
            x = borders_interval[0] if x < borders_interval[0] else x
            x = borders_interval[1] if x > borders_interval[1] else x
            y = round(item[1]  + shift_Y_axis, 1)
            a = -(round(beta, 1))
            f.write(f'n{n} G01 X{x} Y{y} A{a}\n')
            n += 1
    # Откатываемся в стороны для дальнейшей пробивки слоя
    f.write(f'n{n} G01 X{round(thread_points[-1][0] - shift_X_axis - 25.0, 1)} Y{round(thread_points[-1][1] + shift_Y_axis + 150.0, 1)} A{-round(beta, 1) - 360.0}\n')