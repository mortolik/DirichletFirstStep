# 📐 DirichletSolverQt

Интерактивное приложение для численного решения краевой задачи Дирихле для эллиптического уравнения второго порядка методом верхней релаксации (SOR). Реализовано на C++ с использованием Qt 5 и Qt Data Visualization.

---

## 📊 Функциональность

- Решение задачи Дирихле в прямоугольной области методом верхней релаксации
- Поддержка двух задач:
  - ✅ Тестовая задача с точным решением $u^*(x, y) = \sin(\pi x y)$
  - 📈 Основная задача с заданными граничными условиями
- Автоматическое построение:
  - Трёхмерного графика $u(x, y)$
  - Таблиц:
    - $u^*(x, y)$
    - $v(x, y)$ — численное решение
    - $|u^* - v|$ — ошибка
    - $v_2(x, y)$ — уточнённое решение (на сетке $2n \times 2m$)
    - $|v - v_2|$ — оценка погрешности $\varepsilon_2$
- Автоматический подбор оптимального параметра релаксации $\omega$
- Отображение подробной справки:
  - параметры сетки
  - количество итераций
  - невязка на начальном приближении $R(0)$
  - максимальная погрешность $\varepsilon_1$
  - оценка точности $\varepsilon_2$
  - координаты точек максимальной ошибки

---

## 🖥️ Интерфейс

- Вкладка **Тестовая задача** — сравнение численного и точного решения
- Вкладка **Основная задача** — построение решения с оценкой точности по уточнённой сетке
- Визуализация: `Q3DSurface`
- Панель параметров с вводом `n`, `m`, $\varepsilon$, шагов и $\omega$
- Таблицы выводятся отдельно в `QDockWidget`

---

## 📐 Используемые методы

- Разностная схема на прямоугольной сетке
- Метод верхней релаксации (SOR)
- Интерполяционное (смешанное) начальное приближение
- Оценка точности через сравнение с сеткой $2n \times 2m$

---

## 🛠️ Требования

- Qt 5.15+
- Qt Data Visualization module (доступен в Qt Charts/Data tools)
- CMake или `.pro` файл для сборки

---

## 🧪 Запуск

```bash
git clone https://github.com/mortolik/DirichletFirstStep
cd DirichletSolverQt
qmake && make
./DirichletFirstStep
