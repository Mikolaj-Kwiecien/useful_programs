# Rzut ukośny bez oporu i z oporem powietrza

*Raport obliczeniowy: od równań Newtona do numerycznego modelu lotu.*

**Temat:** ruch ciała wyrzuconego pod kątem do poziomu.  
**Metoda:** rozwiązanie analityczne w próżni oraz całkowanie numeryczne dla oporu aerodynamicznego.  
**Dane:** przykład modelowy, nie pomiar laboratoryjny. Wobec braku numeru indeksu przyjęto jawnie $v_0=20\,\mathrm{m/s}$ i $\alpha=\pi/4\,\mathrm{rad}$ (czyli $45^\circ$). Wartości te należy zastąpić dwiema pierwszymi i dwiema ostatnimi cyframi własnego indeksu; kąt podany w stopniach przeliczamy na radiany przed podstawieniem do wzorów.

Rzut ukośny łączy ruch poziomy i pionowy. Bez oporu są to niezależne składowe z jednakowym czasem lotu. Powietrze zmienia prędkość w obu kierunkach, dlatego tor przestaje być parabolą, a pełne rozwiązanie otrzymujemy numerycznie.

## Spis treści

1. [Cel, układ współrzędnych i dane](#1-cel-układ-współrzędnych-i-dane)
2. [Rzut bez oporu powietrza](#2-rzut-bez-oporu-powietrza)
3. [Opór powietrza i równania ruchu](#3-opór-powietrza-i-równania-ruchu)
4. [Wyniki i porównanie modeli](#4-wyniki-i-porównanie-modeli)
5. [Kod i odtwarzanie wykresów](#5-kod-i-odtwarzanie-wykresów)
6. [Założenia i ograniczenia](#6-założenia-i-ograniczenia)
7. [Wnioski](#7-wnioski)

---

## 1. Cel, układ współrzędnych i dane

Początek układu współrzędnych umieszczamy w punkcie wyrzutu. Oś $x$ jest pozioma i skierowana zgodnie z ruchem, a oś $y$ pionowa i skierowana ku górze. Przyspieszenie grawitacyjne ma wartość $g=9{,}81\,\mathrm{m/s^2}$. Ciało startuje z wysokości $y_0=0$ z prędkością początkową $v_0$ pod kątem $\alpha$ do poziomu.

Wektor prędkości początkowej rozkłada się na składowe:

$$
v_{x0}=v_0\cos\alpha,\qquad v_{y0}=v_0\sin\alpha.
$$

| Symbol | Znaczenie | Jednostka SI |
| :--- | :--- | :---: |
| $t$ | czas od wyrzutu | $\mathrm{s}$ |
| $x,y$ | położenie poziome i wysokość | $\mathrm{m}$ |
| $v_x,v_y,v$ | składowe i wartość prędkości | $\mathrm{m/s}$ |
| $a_x,a_y$ | składowe przyspieszenia | $\mathrm{m/s^2}$ |
| $m$ | masa ciała | $\mathrm{kg}$ |
| $g$ | przyspieszenie grawitacyjne | $\mathrm{m/s^2}$ |
| $\rho$ | gęstość powietrza | $\mathrm{kg/m^3}$ |
| $C_D$ | bezwymiarowy współczynnik oporu | $1$ |
| $A$ | pole przekroju ciała | $\mathrm{m^2}$ |
| $\alpha$ | kąt wyrzutu | $\mathrm{rad}$ |
| $R,H,T$ | zasięg, maksymalna wysokość, czas lotu | $\mathrm{m},\mathrm{m},\mathrm{s}$ |

W przykładzie z oporem użyto parametrów zbliżonych do piłki tenisowej: $m=0{,}145\,\mathrm{kg}$, średnica $d=0{,}073\,\mathrm{m}$, $C_D=0{,}47$ oraz $\rho=1{,}225\,\mathrm{kg/m^3}$. Są to przyjęte parametry modelu, a nie uniwersalne stałe dla dowolnego pocisku.

## 2. Rzut bez oporu powietrza

### 2.1. Równania składowych

Jedyną siłą jest ciężar $\vec F_g=(0,-mg)$. Z drugiej zasady Newtona:

$$
m\ddot{x}=0,\qquad m\ddot{y}=-mg,
$$

czyli

$$
\ddot{x}=0,\qquad \ddot{y}=-g.
$$

Całkujemy pierwsze równanie i korzystamy z warunków $v_x(0)=v_0\cos\alpha$, $x(0)=0$:

$$
v_x(t)=v_0\cos\alpha,\qquad x(t)=v_0\cos\alpha\,t.
$$

W kierunku pionowym, przy $v_y(0)=v_0\sin\alpha$ i $y(0)=0$:

$$
v_y(t)=v_0\sin\alpha-gt,
$$

$$
y(t)=v_0\sin\alpha\,t-\frac{1}{2}gt^2.
$$

Po wyeliminowaniu czasu $t=x/(v_0\cos\alpha)$ otrzymujemy równanie toru:

$$
y(x)=x\tan\alpha-\frac{g x^2}{2v_0^2\cos^2\alpha}.
$$

Jest to parabola, ponieważ współczynnik przy $x^2$ jest stały i ujemny.

### 2.2. Wierzchołek toru, czas lotu i zasięg

W najwyższym punkcie $v_y=0$. Z równania prędkości pionowej:

$$
t_H=\frac{v_0\sin\alpha}{g}.
$$

Podstawienie tej chwili do $y(t)$ daje maksymalną wysokość:

$$
H=\frac{v_0^2\sin^2\alpha}{2g}.
$$

Przy lądowaniu na wysokości wyrzutu $y(T)=0$. Nietrywialnym rozwiązaniem równania $v_0\sin\alpha\,T-\tfrac12gT^2=0$ jest:

$$
T=\frac{2v_0\sin\alpha}{g}.
$$

Zasięg to poziome położenie po czasie $T$:

$$
R=x(T)=v_0\cos\alpha\,T
=\frac{v_0^2\sin(2\alpha)}{g}.
$$

Ponieważ dla ustalonego $v_0$ czynnik $v_0^2/g$ jest stały, zasięg maksymalizuje $\sin(2\alpha)$. Jego maksimum wynosi $1$, zatem w tym idealnym modelu $R$ jest największy dla $\alpha=45^\circ$.

### 2.3. Przykład liczbowy

Dla przyjętych $v_0=20\,\mathrm{m/s}$ i $\alpha=\pi/4\,\mathrm{rad}$:

| Wielkość | Obliczenie | Wartość |
| :--- | :--- | ---: |
| $v_{x0}$ | $v_0\cos\alpha$ | $14{,}14\,\mathrm{m/s}$ |
| $v_{y0}$ | $v_0\sin\alpha$ | $14{,}14\,\mathrm{m/s}$ |
| Czas do wierzchołka $t_H$ | $v_{y0}/g$ | $1{,}44\,\mathrm{s}$ |
| Maksymalna wysokość $H$ | $v_{y0}^2/(2g)$ | $10{,}19\,\mathrm{m}$ |
| Czas lotu $T$ | $2v_{y0}/g$ | $2{,}88\,\mathrm{s}$ |
| Zasięg $R$ | $v_{x0}T$ | $40{,}77\,\mathrm{m}$ |

## 3. Opór powietrza i równania ruchu

### 3.1. Model oporu kwadratowego

Dla ciała poruszającego się w powietrzu przyjmujemy siłę aerodynamiczną przeciwną do wektora prędkości, o wartości proporcjonalnej do kwadratu szybkości:

$$
\vec F_D=-\frac{1}{2}\rho C_D A\,|\vec v|\,\vec v,
\qquad |\vec v|=\sqrt{v_x^2+v_y^2}.
$$

Czynnik $|\vec v|\vec v$ zachowuje kierunek siły przeciwny do ruchu i nadaje jej wartość proporcjonalną do $v^2$. Dzieląc równanie Newtona przez masę, definiujemy:

$$
\beta=\frac{\rho C_D A}{2m},\qquad [\beta]=\mathrm{m^{-1}}.
$$

Wówczas układ równań pierwszego rzędu ma postać:

$$
\dot{x}=v_x,\qquad \dot{y}=v_y,
$$

$$
\dot{v}_x=-\beta\sqrt{v_x^2+v_y^2}\,v_x,
$$

$$
\dot{v}_y=-g-\beta\sqrt{v_x^2+v_y^2}\,v_y.
$$

Składowe są sprzężone: opór poziomy zależy również od pionowej składowej szybkości i odwrotnie. Z tego względu wzory na parabolę, zasięg i wysokość z poprzedniej sekcji nie są już rozwiązaniem układu.

### 3.2. Dlaczego energia mechaniczna maleje?

Energię mechaniczną liczymy względem poziomu wyrzutu:

$$
E_m=\frac{1}{2}m(v_x^2+v_y^2)+mgy.
$$

Jej zmiana w czasie jest równa mocy siły oporu. Ciężar jest zachowawczy, a opór wykonuje pracę ujemną:

$$
\frac{dE_m}{dt}=\vec F_D\cdot\vec v
=-\frac{1}{2}\rho C_D A\,|\vec v|^3\leq 0.
$$

Energia mechaniczna układu ciało–Ziemia nie rośnie; część energii przechodzi w energię wewnętrzną powietrza i ciała.

### 3.3. Rozwiązanie numeryczne

Układ integrujemy numerycznie od stanu początkowego $(x,y,v_x,v_y)=(0,0,v_0\cos\alpha,v_0\sin\alpha)$. Integrator `solve_ivp` adaptuje krok obliczeń, a zdarzenie $y=0$ z ruchem w dół zatrzymuje obliczenia przy lądowaniu. Zapisane próbki służą następnie do narysowania wykresów.

## 4. Wyniki i porównanie modeli

Wykresy zapisuje skrypt [generuj_wykresy.py](generuj_wykresy.py) w folderze `files/rzut_ukosny/`. Po jego uruchomieniu poniższe obrazy są osadzane względnymi ścieżkami.

### 4.1. Tor lotu

![Porównanie toru lotu bez oporu i z oporem kwadratowym](files/rzut_ukosny/01_tory_lotu.png)

*Rysunek 1. Opór zmniejsza zasięg i wysokość, a tor nie jest parabolą. Osie podano w metrach.*

### 4.2. Położenie w czasie

![Położenie poziome i wysokość w funkcji czasu](files/rzut_ukosny/02_polozenie_czas.png)

*Rysunek 2. W próżni składowa pozioma rośnie liniowo. Przy oporze prędkość pozioma maleje.*

### 4.3. Składowe prędkości

![Składowe prędkości w funkcji czasu](files/rzut_ukosny/03_predkosc_czas.png)

*Rysunek 3. Bez oporu $v_x$ jest stałe, a $v_y$ zmienia się liniowo. Opór zmienia obie składowe.*

### 4.4. Przyspieszenie

![Składowe przyspieszenia w funkcji czasu](files/rzut_ukosny/04_przyspieszenie_czas.png)

*Rysunek 4. Bez oporu $a_x=0$ i $a_y=-g$. Przy oporze przyspieszenie zależy od chwilowego wektora prędkości.*

### 4.5. Energia mechaniczna

![Energia mechaniczna w funkcji czasu](files/rzut_ukosny/05_energia.png)

*Rysunek 5. W próżni suma energii kinetycznej i potencjalnej jest stała. Opór powoduje jej spadek.*

### 4.6. Zasięg a kąt wyrzutu

![Zasięg w zależności od kąta wyrzutu](files/rzut_ukosny/06_zasieg_kat.png)

*Rysunek 6. W próżni maksimum przypada na $\pi/4\,\mathrm{rad}$. Przy oporze optimum zależy także od parametrów ciała i atmosfery.*

Porównanie ma sens przy tych samych warunkach początkowych. Nie można z góry założyć, że zasięg z oporem ma maksimum dokładnie przy $\pi/4\,\mathrm{rad}$; krzywa w ostatnim wykresie wyznacza je dla przyjętych parametrów.

## 5. Kod i odtwarzanie wykresów

Wykresy przygotowano w Pythonie, korzystając z NumPy, Matplotlib i SciPy. Pełny, uruchamialny skrypt znajduje się w pliku [generuj_wykresy.py](generuj_wykresy.py). Z katalogu głównego repozytorium uruchom:

```powershell
python -m pip install numpy matplotlib scipy
python pages/solutions/04_sol_vscode/generuj_wykresy.py
```

Skrypt zapisuje sześć plików PNG obok raportu, w `pages/solutions/04_sol_vscode/files/rzut_ukosny/`. Parametry wejściowe (`V0`, `ANGLE_DEG`, masa, średnica, $C_D$ i $\rho$) są zebrane na początku programu. Po zmianie dwóch pierwszych wartości trzeba ponownie uruchomić skrypt, aby odświeżyć wykresy. Pełny kod obliczeniowy, wraz z tworzeniem rysunków, jest w pliku `.py`; kluczowy fragment modelu wygląda następująco:

```python
def vacuum_trajectory(speed, angle_deg, times):
	angle = np.deg2rad(angle_deg)
	vx0 = speed * np.cos(angle)
	vy0 = speed * np.sin(angle)
	x = vx0 * times
	y = vy0 * times - 0.5 * G * times**2
	vx = np.full_like(times, vx0)
	vy = vy0 - G * times
	return x, y, vx, vy


def drag_rhs(_time, state):
	_, _, vx, vy = state
	speed = np.hypot(vx, vy)
	return [vx, vy, -BETA * speed * vx, -G - BETA * speed * vy]


def ground_event(_time, state):
	return state[1]


ground_event.terminal = True
ground_event.direction = -1


angle = np.deg2rad(ANGLE_DEG)
solution = solve_ivp(
	drag_rhs,
	(0.0, 30.0),
	[0.0, 0.0, V0 * np.cos(angle), V0 * np.sin(angle)],
	events=ground_event,
	dense_output=True,
	rtol=1e-9,
	atol=1e-11,
	max_step=0.02,
)
```

W praktycznym użyciu skrypt sprawdza powodzenie całkowania, wyznacza próbki aż do zdarzenia lądowania, tworzy wykresy i wypisuje czas lotu, zasięg oraz wysokość maksymalną dla obu modeli.

## 6. Założenia i ograniczenia

- Przyspieszenie grawitacyjne jest stałe, a powierzchnia Ziemi płaska.
- Wysokość wyrzutu i lądowania jest taka sama; nie uwzględniamy przeszkód ani odbicia.
- Powietrze jest nieruchome i ma stałą gęstość $\rho$.
- Opór jest kwadratowy, a współczynnik $C_D$ nie zależy od szybkości ani orientacji ciała.
- Ciało opisano przez stałą masę i stałe pole przekroju; pominięto wiatr, rotację oraz siłę Magnusa.
- Wyniki z oporem są wynikiem numerycznego przybliżenia, a nie zamkniętego wzoru analitycznego.

## 7. Wnioski

1. Bez oporu pozioma składowa prędkości jest stała, pionowa zmienia się pod wpływem stałego przyspieszenia $-g$, a tor jest parabolą.
2. Przy równych wysokościach początku i końca czas lotu wynika z warunku $y(T)=0$, a zasięg z podstawienia $T$ do $x(t)$.
3. Siła oporu kwadratowego zależy od całego wektora prędkości, dlatego ruch poziomy i pionowy stają się sprzężone i wymagają całkowania numerycznego.
4. Opór zmniejsza energię mechaniczną i zazwyczaj obniża wysokość oraz zasięg. Kąt maksymalnego zasięgu zależy wówczas od parametrów ciała i powietrza.
5. Wartości w raporcie są przykładowe. Aby spełnić warunek indywidualny, należy w skrypcie i obliczeniach zastąpić $20\,\mathrm{m/s}$ dwiema pierwszymi cyframi indeksu, a wartość kąta dwiema ostatnimi cyframi indeksu (w stopniach), przeliczając ją na radiany.

---

**Nota:** raport rozdziela model idealny od modelu z oporem, podaje wyprowadzenia analityczne i opis procedury numerycznej. Obrazy wykresów są plikami wynikowymi skryptu, a nie ilustracjami pomiarów.

