# meshgen2
meshgen2 это чрезвычайно простой сеточный генератор. Его код очень компактен. Грубо говоря если вам требуется быстро посчитать обтекание цилиндра то начните с данного простого сеточного генератора.
Генератор meshgen2 написан на языке С/С++. Он выдает данные в формате ASCII .dat Tecplot 360 и Вы можете довольно быстро написать считыватель формата ASCII .dat в своей программе.
Конечно для сложных сеток генератор meshgen2 скорее всего неприменим. Если у Вас такая сетка то сообщите автору (kirill7785@mail.ru) он подумает может быть вашу задачу можно разбить с помощью сеточного генератора meshgen2. 

Характерные примеры простых сеток создаваемых генератором meshgen2 представлены ниже:

![alt_text](https://github.com/kirill7785/meshgen2/blob/main/pic/%D0%9A%D0%BE%D0%BB%D0%B5%D0%BD%D0%BE%D0%BE%D0%B1%D1%80%D0%B0%D0%B7%D0%BD%D1%8B%D0%B9%20%D1%81%D0%BC%D0%B5%D1%81%D0%B8%D1%82%D0%B5%D0%BB%D1%8C.png)
![alt_text](https://github.com/kirill7785/meshgen2/blob/main/pic/%D0%9A%D0%BE%D0%BB%D0%B5%D0%BD%D0%BE%D0%BE%D0%B1%D1%80%D0%B0%D0%B7%D0%BD%D1%8B%D0%B9%20%D1%81%D0%BC%D0%B5%D1%81%D0%B8%D1%82%D0%B5%D0%BB%D1%8C%20%D0%B4%D0%B0%D0%BD%D0%BD%D1%8B%D0%B5.png)
Коленообразный смеситель
![alt_text](https://github.com/kirill7785/meshgen2/blob/main/pic/%D0%9E%D0%B1%D1%82%D0%B5%D0%BA%D0%B0%D0%BD%D0%B8%D0%B5%20%D1%86%D0%B8%D0%BB%D0%B8%D0%BD%D0%B4%D1%80%D0%B0%20.png)
![alt_text](https://github.com/kirill7785/meshgen2/blob/main/pic/%D0%9E%D0%B1%D1%82%D0%B5%D0%BA%D0%B0%D0%BD%D0%B8%D0%B5%20%D1%86%D0%B8%D0%BB%D0%B8%D0%BD%D0%B4%D1%80%D0%B0%20%D0%B4%D0%B0%D0%BD%D0%BD%D1%8B%D0%B5.png)
Обтекание цилиндра.
![alt_text](https://github.com/kirill7785/meshgen2/blob/main/pic/%D0%9A%D0%BB%D0%B5%D0%BD%D0%BE%D0%B2%D1%8B%D0%B9%20%D0%BB%D0%B8%D1%81%D1%82%20.png)
![alt_text](https://github.com/kirill7785/meshgen2/blob/main/pic/%D0%97%D0%B2%D0%B5%D0%B7%D0%B4%D0%B0%20%D0%B4%D0%B0%D0%BD%D0%BD%D1%8B%D0%B5.png)
Кленовый лист.
![alt_text](https://github.com/kirill7785/meshgen2/blob/main/pic/%D0%95%D0%B4%D0%B8%D0%BD%D0%B8%D1%87%D0%BD%D0%B0%D1%8F%20%D0%BE%D0%BA%D1%80%D1%83%D0%B6%D0%BD%D0%BE%D1%81%D1%82%D1%8C%20.png)
![alt_text](https://github.com/kirill7785/meshgen2/blob/main/pic/%D0%9E%D0%BA%D1%80%D1%83%D0%B6%D0%BD%D0%BE%D1%81%D1%82%D1%8C%20%D0%B4%D0%B0%D0%BD%D0%BD%D1%8B%D0%B5.png)
Единичная окружность.
