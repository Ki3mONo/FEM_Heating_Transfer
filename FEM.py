import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.widgets import Slider
import math as Math

##############################################
# Funkcje pomocnicze do wyświetlania wyników #
##############################################

def print_solution(points, solutions):
    print("Rozwiązanie MES:")
    for x_val, u_val in zip(points, solutions):
        print(f"x = {x_val:.6f},  u(x) = {u_val:.6f}")   
        
def show_plot(points,solutions):
    plt.figure(figsize=(16, 9))
    plt.plot(points, solutions)
    plt.title("Równanie transportu ciepła - MES\n Wykres dla funkcji u(x)")
    plt.xlabel("x")
    plt.ylabel("u(x)")
    plt.grid(True)
    plt.gca().get_lines()[0].set_color("green")
    plt.show()

def save_solution_to_file(points, solutions, filename):
    with open(filename, "w") as file:
        file.write("x,u(x)\n")
        for x_val, u_val in zip(points, solutions):
            file.write(f"{x_val},{u_val}\n")

def show_plot_with_table(filename="solution.csv", rows_per_view=50):

    df = pd.read_csv(filename).round(6)
    points = df["x"]
    solutions = df["u(x)"]

    fig = plt.figure(figsize=(16, 9))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.80, 0.19, 0.01])

    ax_plot   = fig.add_subplot(gs[0, 0])
    ax_table  = fig.add_subplot(gs[0, 1])
    ax_slider = fig.add_subplot(gs[0, 2])
    
    ax_plot.plot(points, solutions)
    ax_plot.set_title("Równanie transportu ciepła - MES\n Wykres i tabela wyników dla funkcji u(x)")
    ax_plot.set_xlabel("x")
    ax_plot.set_ylabel("u(x)")
    ax_plot.grid(True)
    
    ax_table.axis("off")
    def update_table(start_idx):
        ax_table.clear()
        ax_table.axis("off")

        end_idx = min(start_idx + rows_per_view, len(df))
        subset = df.iloc[start_idx:end_idx]

        tbl = ax_table.table(
            cellText=subset.values,
            colLabels=subset.columns,
            loc="center"
        )
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(8)

        plt.draw()
        
    update_table(0)

    max_slider = max(0, len(df) - rows_per_view)

    slider = Slider(
        ax=ax_slider,
        label="",
        orientation='vertical', 
        valmin=0,                
        valmax=max_slider,       
        valinit=max_slider,     
        valstep=1,            
    )
    slider.poly.set_visible(False)
    slider.valtext.set_visible(False)

    def on_slider_change(val):
        start_idx = max_slider - int(val)
        update_table(start_idx)

    slider.on_changed(on_slider_change)

    ax_plot.get_lines()[0].set_color("green")
    plt.tight_layout()
    plt.show()


##################################################################
# Funkcje do całkowania dwupunktową kwadraturą Gaussa-Legendre'a #
##################################################################

def gauss_legendre_quadrature(func, start, end, points, weights):
    """
    Całkuje funkcję func(x) na przedziale [start,end] 
    dwupunktową kwadraturą Gaussa-Legendre'a.
    """
    # Obliczamy jakobian przejścia i ustalamy wartość początkową całki
    # Domyślnie kwadratura Gaussa-Legendre'a jest zdefiniowana na przedziale [-1,1]
    jacobian = 0.5 * (end - start)
    integral = 0.0
    
    # Obliczamy wartość całki przez sumowanie wartości w punktach Gaussa
    for point, weight in zip(points, weights):
        # Mapujemy punkt Gaussa na przedział [start,end]
        x_mapped = 0.5*(start+end) + jacobian*point
        # Dodajemy wartość w punkcie Gaussa pomnożoną przez wagę
        integral += weight * func(x_mapped)
    
    # Zwracamy wartość całki pomnożoną przez jakobian
    return integral * jacobian

def get_gauss_legendre_points_and_weights():
    """
    Zwraca punkty i wagi dwupunktowej kwadratury 
    Gaussa-Legendre'a na przedziale referencyjnym [-1,1].
    """
    gauss_points = np.array([-1.0 / Math.sqrt(3), 1.0 / Math.sqrt(3)])
    gauss_weights = np.array([1.0, 1.0])
    return gauss_points, gauss_weights


#####################################
# Część związana z rozwiązaniem MES #
#  dla równania transportu ciepła   #
#####################################

def generate_points(start, end, number_of_elements):
    """
    Generuje równomiernie rozłożone punkty (węzły).
    Przedział [start, end] dzielony jest na number_of_elements części.
    """
    return np.linspace(start, end, number_of_elements + 1)

def e_i(node_positions, element_index):
    """
    Zwraca krotkę (x_left, x_right, element_length) dla
    i-tego elementu skończonego w tablicy node_positions,
    gdzie:
        x_left        - współrzędna lewego węzła,
        x_right       - współrzędna prawego węzła,
        element_length- różnica x_right - x_left.
    """
    x_left = node_positions[element_index]
    x_right = node_positions[element_index + 1]
    element_length = x_right - x_left
    return x_left, x_right, element_length

def generate_B_matrix(node_positions):
    """
    Funkcja tworząca macierz globalną (tzw. macierz B) dla zagadnienia
    z warunkiem brzegowym typu Robina (na lewym węźle) i Dirichleta (na prawym węźle).
    """
    # Pobieramy liczbę węzłów i elementów
    number_of_nodes = len(node_positions)
    number_of_elements = number_of_nodes - 1

    # Inicjalizacja macierzy B
    B_matrix = np.zeros((number_of_nodes, number_of_nodes))

    # Pobieramy punkty i wagi Gaussa
    gauss_points, gauss_weights = get_gauss_legendre_points_and_weights()

    # Petla po elementach skończonych
    for element_index in range(number_of_elements):
        # Pobieramy współrzędne lewego i prawego węzła oraz długość elementu
        left_node, right_node, element_length = e_i(node_positions, element_index)

        # Lokalna macierz sztywności (2x2) bo każdy element ma dwa węzły
        local_matrix = np.zeros((2, 2))

        # Obliczamy lokalną macierz stosując kwadraturę Gaussa
        for local_i in range(2):
            for local_j in range(2):
                # Obliczamy wartość funkcji podcałkowej dla i-tego i j-tego wiersza
                def integrand(x):
                    # Pochodne funkcji kształtu w rzeczywistej dziedzinie elementu
                    shape_function_derivatives_domain = np.array([-1.0 / element_length, 1.0 / element_length])
                    # Wartość funkcji k w punkcie x 
                    k = K(x)
                    # Zwracamy wartość funkcji podcałkowej
                    return (k*shape_function_derivatives_domain[local_i]*shape_function_derivatives_domain[local_j])
                # Obliczamy wartość całki z funkcji podcałkowej
                local_matrix[local_i, local_j] = gauss_legendre_quadrature(
                    integrand, left_node, right_node, gauss_points, gauss_weights
                )

        # Wprowadzamy wartości lokalnej macierzy sztywności do macierzy globalnej
        global_indices = [element_index, element_index + 1]
        for local_i in range(2):
            for local_j in range(2):
                # Dodajemy wartość z lokalnej macierzy do globalnej
                B_matrix[global_indices[local_i], global_indices[local_j]] += local_matrix[local_i, local_j]

        
    # Warunek brzegowy typu Robina na węźle 0 (dodanie 1.0 na elemencie diagonalnym)
    B_matrix[0, 0] += 1.0

    # Warunek brzegowy typu Dirichleta na ostatnim węźle:
    # "Zerujemy" cały wiersz i ustawiamy 1.0 na elemencie diagonalnym
    B_matrix[(number_of_nodes - 1), :] = 0.0
    B_matrix[(number_of_nodes - 1), (number_of_nodes - 1)] = 1.0

    # Zwracamy macierz B
    return B_matrix

def generate_L_matrix(node_positions):
    """
    Funkcja tworząca globalny wektor prawej strony dla zagadnienia
    z warunkiem brzegowym typu Robina (na lewym węźle) i
    Dirichleta (na prawym węźle).
    """
    # Pobieramy liczbę węzłów i elementów
    number_of_nodes = len(node_positions)
    number_of_elements = number_of_nodes - 1

    # Inicjalizacja macierzy L
    L_vector = np.zeros(number_of_nodes)

    # Pobieramy punkty i wagi Gaussa
    gauss_points, gauss_weights = get_gauss_legendre_points_and_weights()

    # Pętla po elementach skończonych
    for element_index in range(number_of_elements):
        # Pobieramy współrzędne lewego i prawego węzła oraz długość elementu
        left_node, right_node, element_length = e_i(node_positions, element_index)

        # Lokalny wektor obciążenia, 2-elementowy bo każdy element z dwiema funkcjami kształtu ma wkłady do 2 węzłów
        local_load_vector = np.zeros(2)

        # Obliczamy lokalny wektor obciążenia stosując kwadraturę Gaussa
        for local_i in range(2):
            # Definiujemy funkcję podcałkową odpowiadającą wkładowi do lokalnego wektora obciążenia
            def load_integrand(x):
                # Wartości funkcji kształtu w rzeczywistej dziedzinie elementu
                shape_functions = np.array([
                    1.0 - (x - left_node) / element_length,  # Funkcja kształtu przypisana do pierwszego węzła
                    (x - left_node) / element_length         # Funkcja kształtu przypisana do drugiego węzła
                ])
                # Obliczamy wartość funkcji podcałkowej jako iloczyn funkcji wymuszenia f(x) i funkcji kształtu
                return f(x) * shape_functions[local_i]

            # Całkujemy funkcję podcałkową metodą Gaussa-Legendre’a w celu wyznaczenia wkładu do lokalnego wektora obciążenia
            local_load_vector[local_i] = gauss_legendre_quadrature(
                load_integrand, left_node, right_node, gauss_points, gauss_weights
            )

        # Wprowadzamy wartości lokalnego wektora obciążenia do globalnego wektora obciążeń
        global_indices = [element_index, element_index + 1]
        for local_i in range(2):
            L_vector[global_indices[local_i]] += local_load_vector[local_i]

    # Warunek brzegowy typu Robina na węźle 0:
    # Odejmujemy wartość warunku brzegowego od pierwszego elementu wektora L
    L_vector[0] -= ROBIN_START_CONDITION

    # Warunek Dirichleta na ostatnim węźle:
    # Ustawiamy wartość warunku brzegowego na ostatnim elemencie wektora L
    L_vector[number_of_nodes - 1] = DIRICHLET_END_CONDITION

    # Zwracamy wektor L
    return L_vector

def solve_heat_equation(start, end, number_of_elements):
    """
    Rozwiązuje równanie transportu ciepła:
       -d/dx [k(x) d/dx u(x)] = 100*x^2
    na [start, end], z użyciem number_of_elements elementów skończonych.
    z warunkiem:
       - Robina w x=0 -> u'(0) + u(0) = 20
       - Dirichleta w x=end -> u(end) = -20
    """
    # Generujemy siatkę węzłów
    node_positions = generate_points(start, end, number_of_elements)
    # Składamy układ równań MES z użyciem siatki węzłów
    B_matrix = generate_B_matrix(node_positions)
    L_matrix = generate_L_matrix(node_positions)
    # Rozwiązujemy układ równań
    solution = np.linalg.solve(B_matrix, L_matrix)
    # Zwracamy punkty i rozwiązanie
    return node_positions, solution

#######################################################
# Ustawienie warunków brzegowych i funkcji k(x), f(x) #
#               u'(0) + u(0) = 20                     #
#                  u(2) = -20                         #   
#######################################################

ROBIN_START_CONDITION = 20.0
DIRICHLET_END_CONDITION = -20.0

def K(x):
    """
    Funkcja k(x) = 1 dla x <= 1, 
                   2*x dla x > 1.
    """
    if x <= 1.0:
        return 1.0
    else:
        return 2.0 * x

def f(x):
    """
    Funkcja f(x) = 100*x^2.
    """
    return 100 * x**2.

#########################
# Główna część programu #
#########################
if __name__ == "__main__":
    
    # Pobieramy liczbę elementów skończonych od użytkownika
    number_of_elements = int(input("Podaj liczbę elementów skończonych: "))

    # Usatwiamy przedział [0,2] jako dziedzinę 
    start_domain = 0.0
    end_domain   = 2.0
    
    # Rozwiązujemy zadane równanie transportu ciepła
    points, solutions = solve_heat_equation(start_domain, end_domain, number_of_elements)

    # Uproszczone wyświetlenie wyników
    print_solution(points, solutions)
    show_plot(points, solutions)
    
    # Rozszerzone Wyświeetlenie wyników ze sliderem tabelki (x,u(x))
    #   Aby z niej skorzystać, zakomentuj powyższe linie i odkomentuj poniższe
    #   Uwaga: zapisze plik "solution.csv" z wynikami w bieżącym katalogu
    
    # save_solution_to_file(points, solutions, "solution.csv")
    # show_plot_with_table("solution.csv")

