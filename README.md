# Repository 1
 import matplotlib.pyplot as plt

class Physic:
    def __init__(self, Z, A):
      '''
        Инициализация и рассчёт основных значений :
        :param Z:  Type is Integer. Зарядовое число элемента.
        :param A: Type is Integer. Массовое число элемента.
        :param N: Type is Integer. Число нейтронов в ядре.
        :param R: Type is float. Радиус ядра.
        :param E: Type is float. Энергия связи, рассчитанная по формуле Вайзекера
        :param e: Type is float. Удельная энергия связи.
        :param M Type is float. Масса атома.
        '''
      self.Z = Z
      self.A = A
      self.N=self.Z-self.A
      self.R=1.2*self.A**(1/3) #в Фм=10**-13 см
      if (self.Z%2==0 and self.A%2==0 ):
        self.E=15.75*self.A - 17.8*self.A**(2/3) - 0.711*self.Z**2*self.A**(-1/3) - 23.7*(self.A/2 - self.Z)**2/self.A + 34*self.A**(-3/4)
      elif (self.Z%2!=0 and self.A%2!=0 ):
        self.E=15.75*self.A - 17.8*self.A**(2/3) - 0.711*self.Z**2*self.A**(-1/3) - 23.7*(self.A/2 - self.Z)**2/self.A - 34*self.A**(-3/4)
      else:
        self.E=15.75*self.A - 17.8*self.A**(2/3) - 0.711*self.Z**2*self.A**(-1/3) - 23.7*(self.A/2 - self.Z)**2/self.A
      self.M=self.E-self.Z*938.272 -self.N*939.565 #масса в МэВ
      self.e=self.E/self.A # удельная энергия связи

    def radius(self, A):
      '''
        Устанавливает значение радиуса для заданного массового числа . Вынесена в отдельную функцию для удобства при построении графика 1.
        :param A: Type is Integer. Массовое число элемента.
        :return: Type is float. Радиус ядра для заданного массового числа.
        '''
      return 1.2*A**(1/3)

    def energy(self, Z):
      '''
        Устанавливает значение энергии для устойчивых ядер (A=2*Z) с зарядовым числом от 1 до 20. Вынесена в отдельную функцию для удобства при построении графика 2.
        Так как А=2*Z возможны только чётно-чётное и нечётное состояния (пятый член в формуле положительный или равен 0, соответственно).
        Четвёртый член зануляется т.к. содержит множитель (А-2*Z).
        :param Z:  Type is Integer. Зарядовое число элемента.
        :param A: Type is Integer. Массовое число элемента.
        :return: Type is float. Энергия связи для заданного зарядового числа
        '''
      A=2*Z
      if (Z%2==0 ):
        return 15.75*A - 17.8*A**(2/3) - 0.711*Z**2*A**(-1/3) + 34*A**(-3/4)
      else:
       return 15.75*A - 17.8*A**(2/3) - 0.711*Z**2*A**(-1/3)


    def beta_decay(self):
      '''
       Проверяет устойчив ли изотоп к бета распаду. Выводит результат проверки в текстовом формате на экран. Проверяет возможно ли деление изотопа на два чётно-чётных осколка. Выводит результат проверки в текстовом формате на экран.
       :return: ничего не возращается
        '''
      element_у = Physic(self.Z+1, self.A);
      if (self.M > element_у.M):
        print ('не устойчив к бета распаду')
      else:
         print ('устойчив к бета распаду')

    def twins_decay(self):
      '''
       Проверяет возможно ли деление изотопа на два чётно-чётных осколка. Выводит результат проверки в текстовом формате на экран.
       :return: ничего не возращается
        '''
      if ( self.Z % 4 == 0 and self.A % 4 == 0):
        print('Возможно деление на 2 одинаковых четно-четных осколка')
      else:
        print('Не возможно деление на 2 одинаковых четно-четных осколка')

    def plot1(self):
      '''
        Функция строит график зависимости радиуса атома от массового числа
        :return: Ничего не возвращается
        '''
      arA = list(range(1,261))
      arR = list()
      for A in arA:
        arR.append(self.radius(A))

      plt.figure(figsize=[9,6])
      plt.plot(arA, arR, linewidth=2)
      plt.xlabel('Массовое число (А)')
      plt.ylabel('Радиус атома')
      plt.title('Зависимость радиуса атома от массового числа')
      plt.grid(True)
      plt.show()

    def plot2(self):
      '''
        Функция строит график зависимости энергии атома от атомного номера для устойчивого состояния первых 20-ти элементов таблицы Менделеева
        :return: Ничего не возвращается
        '''
      arZ = list(range(1,21))
      arE = list()
      for Z in arZ:
        arE.append(self.energy(Z))

      plt.figure(figsize=[4,6])
      plt.plot(arZ, arE, linewidth=2)
      plt.xlabel('Зарядовое число (Z)')
      plt.ylabel('Энергия устойчивого атома')
      plt.title('Зависимость энергии атома от атомного номера')
      plt.grid(True)
      plt.show()
U=Physic(92,238)
print("удельная энергия связи =", U.E, "; масса атома =", U.M, "; радиус атома =", U.R)
U.beta_decay()
U.twins_decay()
U.plot1()
U.plot2()

