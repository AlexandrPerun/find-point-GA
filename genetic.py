import math
import random
import time
from haversine import haversine

def _get_azimuth(point_1, point_2):
    '''
    Method for determining the azimuth - the angle between the north direction and the direction from the start point
    to the endpoint
    :param point_1: coordinates of the starting point (latitude, longitude)
    :param point_2: coordinates of the endpoint (latitude, longitude)
    :return: azimuth in degrees
    '''

    # in radians
    lat1 = math.radians(point_1[0])
    lat2 = math.radians(point_2[0])
    long1 = math.radians(point_1[1])
    long2 = math.radians(point_2[1])
    # cosines and sines of latitudes and difference of longitudes
    cos_lat1 = math.cos(lat1)
    cos_lat2 = math.cos(lat2)
    sin_lat1 = math.sin(lat1)
    sin_lat2 = math.sin(lat2)

    delta = long2 - long1
    cdelta = math.cos(delta)
    sdelta = math.sin(delta)

    # calculation of the initial azimuth
    x = (cos_lat1 * sin_lat2) - (sin_lat1 * cos_lat2 * cdelta)
    y = sdelta * cos_lat2
    z = math.degrees(math.atan(-y / x))

    if (x < 0):
        z = z + 180.

    z2 = (z + 180.) % 360. - 180.

    z2 = - math.radians(z2)
    anglerad2 = z2 - (2 * math.pi) * math.floor(z2 / (2 * math.pi))
    angledeg = math.degrees(anglerad2)
    return  angledeg

# входные данные
B=[46.1578860816622, -1.13280177915811] # точка - начало известного вектора
A=[46.1574781369847, -1.13255540981461] # точка - конец известного вектора
distance = 46044.03676089714 # длина неизвестного вектора в миллиметрахб
angle = -89.43474295596019 # угол между неизвестныи и известным векторами (-180, 180)

azimuth_1 = _get_azimuth(B, A) # азимут известного вектора
if azimuth_1 > 180:
    azimuth_1 -= 360

azimuth2 = azimuth_1 + angle # вичисление азимута неизвестного вектора
if azimuth2 < -180:
    azimuth2 += 360
elif azimuth2 > 180:
    azimuth2 -= 360

quantity_of_chrom = 100 # колличество особей в популяции
probability_of_mutation = 0.1 # вероятность мутации

# область поиска неизвестной точки С (приращение точки В на +- 100 километров
long_min = (B[1] - 0.001)
long_max = (B[1] + 0.001)
lat_min = (B[0] - 0.001)
lat_max = (B[0] + 0.001)





#количество битов для каждой переменной
def exponent_m(x_min, x_max):
    m = (x_max - x_min)*10**9 #10**5 так как выбранная нами точность 5 знаков после запятой
    m = math.log(m, 2)
    m = math.ceil(m)
    return m

#генератор начальной популяции
def chromosomeGenerator(m1, m2, quantity_of_chrom):
    random_gen1 = ''
    random_gen2 = ''
    population = []
    for i in range(quantity_of_chrom):
        for j in range(m1):
            gen = random.randint(0, 1)
            gen = str(gen)
            random_gen1 += gen

        for j in range(m2):
            gen = random.randint(0, 1)
            gen = str(gen)
            random_gen2 += gen

        chromosome = random_gen1 + random_gen2
        population.append(chromosome)
        random_gen1 = ''
        random_gen2 = ''
    return population

#декодер популяции из генотипа в фенотип
def decoder(population, lat_min, lat_max, long_min, long_max, m1, m2):
    population_dec = []
    for i in range(len(population)):
        chromosome = population[i]
        bin_lat = chromosome[:m1]
        bin_long = chromosome[m1:]
        dec_lat = int(bin_lat, 2)
        dec_long = int(bin_long, 2)
        lat = lat_min + dec_lat * (lat_max - lat_min) / (2 ** m1 - 1)
        long = long_min + dec_long * (long_max - long_min) / (2 ** m2 - 1)
        population_dec.append([lat, long])
    return population_dec

# целевая функция
def eval_function(population_dec, distance, B, azimuth2):
    eval_list = []
    for i in range(len(population_dec)):
        lat = population_dec[i][0]
        long = population_dec[i][1]
        C = [lat, long]
        angle = _get_azimuth(B, C)
        k = abs(azimuth2 - angle)/10000
        g = abs(distance - (haversine(B, C) * 1000000))/1000000
        f = -(g+k)
        eval_list.append(f)
    return eval_list

#вычисление совокупной вероятности
def cumulative_probability(eval_list):
    f = 0
    probability_list = []
    q = 0
    cumulative_probability_list = []
    for i in eval_list:
        f += i - min(eval_list) #общая функция соответствия
    for i in eval_list:
        if f != 0:
            p = ((i - min(eval_list))/f)    #вероятность отбора каждой хромосомы
        else:
            p = 1 / len(eval_list)
        probability_list.append(p)
    for i in probability_list:
        q += i #совокупная вероятность
        cumulative_probability_list.append(q)
    return cumulative_probability_list

#"колесо рулетки" для отбора новой популяции
def roulette_wheel(population, q_list):
    new_population = []
    for i in range(len(population)):
        r = random.random()
        for j in range(len(population)):
            if r <= q_list[j]:
                new_population.append(population[j])
                break
    return new_population

# скрещивание
def crossbreeding(population):
    cross_population = []
    index_list = []
    rand_position = random.randint(1, len(population[0])-1)
    new_cross_population = []
    # генерируем списки со случайными хромосомами и их индексами
    for i in range(len(population)):
        r = random.random()
        if r < 0.25:  # вероятность скрещивания
            cross_population.append(population[i])
            index_list.append(i)
    # скрещивание выбраных хромосом
    m = len(cross_population)
    if m > 1:
        for i in range(0, len(cross_population) - 1, 2):
            parent1 = cross_population[i]
            parent2 = cross_population[i + 1]
            chrom1 = parent1[:rand_position] + parent2[rand_position:]
            chrom2 = parent2[:rand_position] + parent1[rand_position:]
            new_cross_population.append(chrom1)
            new_cross_population.append(chrom2)
        new_cross_population.append(cross_population[len(cross_population) - 1])
        # замена хромосом на скрещенные
        j = 0
        for i in index_list:
            population[i] = new_cross_population[j]
            j += 1
    return population

#мутация генов
def mutation(population, probability_of_mutation):
    chrom_len = len(population[0])
    quantity_gen = chrom_len * len(population)
    index_gen = []
    for i in range(1, quantity_gen):
        r = random.random()
        if r < probability_of_mutation:
            index_gen.append(i)
    for i in index_gen:
        number_chrom = i // chrom_len
        gen_in_chrom = i % chrom_len - 1
        if gen_in_chrom == -1:
            number_chrom -= 1
            gen_in_chrom = chrom_len - 1

        chromosome = population[number_chrom]
        gen = chromosome[gen_in_chrom]
        if gen == '0':
            gen = '1'
        else:
            gen = '0'
        mut_chrom = chromosome[:gen_in_chrom] + gen + chromosome[gen_in_chrom + 1:]
        population[number_chrom] = mut_chrom

    return population

def prog():
    start_time = time.time()
    m1 = exponent_m(lat_min, lat_max)
    m2 = exponent_m(long_min, long_max)
    first_population = chromosomeGenerator(m1, m2, quantity_of_chrom)
    population = list(first_population)
    best_chrom_list =[]
    n = 300 #максимальное кол-во итераций для вывода наилучшей особи
    for i in range(n):
        population_dec = decoder(population, lat_min, lat_max, long_min, long_max, m1, m2)
        eval_list = eval_function(population_dec, distance, B, azimuth2)

        index_best_chrom = eval_list.index(max(eval_list))
        best_chrom = population_dec[index_best_chrom]
        best_chrom_list.append([max(eval_list), best_chrom, i+1])

        q_list = cumulative_probability(eval_list)
        select_pop = roulette_wheel(population, q_list)
        population = list(select_pop)
        cross_pop = crossbreeding(population)
        population = list(cross_pop)
        mut_pop = mutation(population, probability_of_mutation)
        population = list(mut_pop)

    time_of_work = (time.time() - start_time)*1000
    time_of_work = round(time_of_work, 5)
    best_of_the_best = max(best_chrom_list)
    print('\nЛучшая особь - {0}, была получена в {1} поколении. Значение целевой функции - {2}'.format(best_of_the_best[1], best_of_the_best[2], best_of_the_best[0]))
    print('Количество особей в поколении: {0}; вероятность мутации: {1}'.format(quantity_of_chrom, probability_of_mutation))
    print('Время работы программы: {0} миллисекунд'.format(time_of_work))

prog()