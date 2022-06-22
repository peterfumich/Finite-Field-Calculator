import tkinter as tk
from matplotlib import pyplot as plt
import numpy as np
from PIL import Image, ImageTk
parameters = {'field_size':4,'modulus':2}
keys = ['0']
values = ['0']
for i in range(parameters['field_size']):
    keys.append(str(i + 1))
    values.append('a^' + str(i))

dictionary = dict(zip(keys, values))

reverse_dictionary = dict(zip(values, keys))
dictionaries = {'dictionary':dictionary,'reverse_dictionary':reverse_dictionary,'values':values}
# IMPORT DATA
def poly_array(poly_short_array):
    length = len(poly_short_array)
    out_array = []
    for coef in poly_short_array:
        new_coef = np.zeros(length)
        new_coef[0] = coef
        out_array.append(new_coef)
    return(out_array)
def Build_Conway_Table():
    with open('CPimport.txt') as file:
        lines = file.readlines()
    coef_row = []
    lookup_row = []
    current_coef = 2
    global coefficients_out
    coefficients_out = []
    global lookup
    lookup = []
    for line in lines:
        if line[0] == "[":
            entry = line[0:-2]
            #Begin parsing new value into table:

            array_entry = eval(entry)
            parameters = [array_entry[0],array_entry[1]]
            coefficients = poly_array(array_entry[2])
            if parameters[0]!=current_coef:
                current_coef = parameters[0]
                lookup.append(lookup_row)
                coefficients_out.append(coef_row)

                lookup_row = []
                coef_row = []

                coef_row.append(coefficients)
                lookup_row.append(parameters)
            else:
                coef_row.append(coefficients)
                lookup_row.append(parameters)
Build_Conway_Table()
def Conway_Polynomial(modulus, power):
    index = lookup[modulus-2].index([modulus, power])
    return(coefficients_out[modulus-2][index])
def Update_Dicts(new_field_size):
    parameters['field_size'] = new_field_size
    keys = ['0']
    values = ['0']
    for i in range(new_field_size):
        keys.append(str(i + 1))
        values.append('a^' + str(i))
    dictionaries['dictionary'] = dict(zip(keys, values))
    dictionaries['reverse_dictionary'] = dict(zip(values, keys))
    dictionaries['values'] = values
def simplify_key(coef_key):
    field_size = parameters['field_size']
    modulus = parameters['modulus']
    depth = int(2 * field_size)
    power = np.log(field_size)/np.log(modulus)
    conway_poly_root = Conway_Polynomial(modulus, power)
    conway_poly = Poly(conway_poly_root)
    print('conway=')
    conway_poly.Get_polystr()
    conway_key = [coef[0] for coef in conway_poly.keyed_coef]
    original_coef_key = coef_key
    if np.sum(np.array(original_coef_key)) <= 1:
        return (original_coef_key)
    shifted_summed_conway_keys = [[[conway_key[(k - (j + 1))] for k in range(field_size)] for j in range(field_size)]]
    for i in range(depth):
        new_layer = []
        for ii in range(field_size ** (i + 1)):
            last_value = shifted_summed_conway_keys[i][ii]
            for j in range(field_size):
                shifted_conway_key = [conway_key[(k - (j + 1))] for k in range(field_size)]
                new_value = [(last_value[k] + shifted_conway_key[k]) % modulus for k in range(field_size)]
                if np.sum((np.array(new_value) + original_coef_key) % modulus) != 1:
                    new_layer.append(new_value)
                else:
                    return((np.array(new_value) + original_coef_key) % modulus)
        shifted_summed_conway_keys.append(new_layer)
def String_to_Poly(string):
    field_size = parameters['field_size']
    local_string = ''
    generate_coef_array = False
    coef_array = np.zeros(field_size)
    out_array = np.zeros((field_size, field_size))
    reverse_dictionary = dictionaries['reverse_dictionary']
    for char in string:
        if char == "(":
            generate_coef_array = True
            coef_array = np.zeros(field_size)
            local_string = ''
        elif char == ")":
            coef_array[int(reverse_dictionary[local_string]) - 1] = 1
            local_string = ''
            generate_coef_array = False
        elif char == "+":
            if generate_coef_array == True:
                coef_array[int(reverse_dictionary[local_string]) - 1] = 1
            else:
                if local_string != '':
                    power = int(local_string[-1])
                else:
                    power = 0
                out_array[power] = coef_array
            local_string = ''
        else:
            local_string += char
    if local_string != '':
        power = int(local_string[-1])
    else:
        power = 0
    out_array[power] = coef_array
    return (out_array)
class Poly():
    def __init__(self, in_keyed_coef):
        field_size = parameters['field_size']
        self.eval_points = []
        self.polystr = ''
        self.keyed_coef = []

        self.coef = np.empty(field_size, dtype=object)

        if type(in_keyed_coef) == str:
            for i in range(field_size):
                self.keyed_coef.append(np.random.randint(0, 2, size=field_size))
        else:
            while len(in_keyed_coef)<field_size:
                in_keyed_coef.append(np.zeros(field_size))
            for coef_key in in_keyed_coef:
                coef_key = coef_key.tolist()
                while len(coef_key)<field_size:
                    coef_key.append(0)
            self.keyed_coef = in_keyed_coef
        self.Set_Coef()
        self.Set_polystr()

    def Rebuild(self,new_parameters):#parameters = [new_modulus,new_power]
        field_size = parameters['field_size']
        new_modulus, new_field_size = new_parameters
        local_f = self.keyed_coef
        new_f = np.zeros((new_field_size, new_field_size))
        parameters['modulus'] = new_modulus
        if field_size <= new_field_size:
            for i in range(field_size):
                for j in range(field_size):
                    new_f[i][j] += local_f[i][j]
        else:
            for n, coef_array in enumerate(local_f):
                for m, coef in enumerate(coef_array):
                    new_f[n % new_field_size][m % new_field_size] += coef
                    new_f[n % new_field_size][m % new_field_size] = new_f[n % new_field_size][m % new_field_size] % new_modulus
        Update_Dicts(new_field_size)
        polynomials['f'] = Poly(new_f)
        self.__init__(polynomials['f'].Get_Keyed_Coef())

    def Set_Coef(self):
        dictionary = dictionaries['dictionary']
        for n, coef_key in enumerate(self.keyed_coef):
            coef = []
            for m, j in enumerate(coef_key):
                if j == 1:
                    coef.append(dictionary[str(m + 1)])
                else:
                    coef.append('0')
            self.coef[n] = coef

    def Set_polystr(self):
        string = ''
        for power, coefficient in enumerate(self.coef):
            string += '+('
            for j in coefficient:
                if j != '0':
                    string += (j + "+")

            if string[-2:] != '+(':
                string = string[:-1]
                string += ')'
                if power != 0:
                    string += 'x^'
                    string += str(power)
            else:
                string = string[:-2]
        self.polystr = string[1:]

    def Get_polystr(self):
        print('Polynomial :', str(self.polystr))
        return (self.polystr)

    def Get_Keyed_Coef(self):
        print('Keyed Values: ', self.keyed_coef)
        return (self.keyed_coef)

    def Get_Coef(self):
        print('Coef :', self.coef)
        return (self.coef)

    def Eval_Points(self,points):
        field_size = parameters['field_size']
        print('evaluating points')
        eval_points = []
        values = dictionaries['values']
        for k, x in enumerate(points):
            fx = np.zeros(field_size)
            for n,coef in enumerate(self.keyed_coef):
                for m,key in enumerate(coef):
                    if key != 0:
                        if x !='0':
                            fx[(k*n+m)%field_size] = 1
                        elif n==0:
                            fx[m] = 1
                        break
            fx = simplify_key(fx)
            new_point = 0
            for kk, key in enumerate(fx):
                if key != 0:
                    new_point = (kk+1)%(parameters['field_size']+1)
            #eval_points.append(values[new_point])
            eval_points.append(new_point)

        return(eval_points)
    def Add(self, poly):  # Adds poly to self(rewriting coef for self)#NOT CORRECT

        sum = (np.array(self.keyed_coef) + np.array(poly.keyed_coef)) % parameters['modulus']
        return Poly(sum)

    def Multiply(self, poly):
        modulus = parameters['modulus']
        field_size = parameters['field_size']
        keyed_coef = np.array(self.keyed_coef)
        keyed_coef[...] = 0
        for n, coef_array in enumerate(self.keyed_coef):
            for m, poly_key_coef in enumerate(poly.keyed_coef):
                key = np.zeros(field_size)
                for a_index, a_coef_key in enumerate(coef_array):
                    for b_index, b_coef_key in enumerate(poly_key_coef):
                        key[(a_index + b_index) % field_size] += (a_coef_key * b_coef_key)
                        key[(a_index + b_index) % field_size] = int(key[(a_index + b_index) % field_size] % modulus)
                for k, key_value in enumerate(key):
                    keyed_coef[(n + m) % field_size][k] += key_value
                    keyed_coef[(n + m) % field_size][k] = keyed_coef[(n + m) % field_size][k] % modulus
        return (Poly(keyed_coef))

    def Simplify(self):

        field_size = parameters['field_size']
        #First get the appropriate Conway Polynomial
        power = np.log(field_size)/np.log(parameters['modulus'])
        if power-int(power)==0:
            new_keys = []
            for coef_key in self.keyed_coef:
                coef_key = simplify_key(coef_key)
                new_keys.append(coef_key)
            self.keyed_coef = new_keys
            self.Set_Coef()
            self.Set_polystr()
        else:
            print('invalid')
polynomials = {'f':Poly, 'g': Poly, 'h': Poly}
def close():
    sand_window.quit()


def Randomize():
    parameters['field_size'] = field_size_tk.get()
    Update_Dicts(parameters['field_size'])
    new_f = Poly('random')
    f.set(new_f.Get_polystr())
    new_g = Poly('random')
    g.set(new_g.Get_polystr())
    polynomials['f'] = new_f
    polynomials['g'] = new_g


def update_function_f(*args):
    polynomials['f'] = Poly(String_to_Poly(f_entry.get()))


def update_function_g(*args):
    polynomials['g'] = Poly(String_to_Poly(g_entry.get()))


def update_functions(*args):
    modulus = parameters['modulus']
    new_field_size = int(field_size_tk.get())
    Update_Dicts(new_field_size)
    local_f = polynomials['f']
    local_f.Rebuild([modulus,new_field_size])
    local_g = polynomials['g']
    local_g.Rebuild([modulus,new_field_size])
    polynomials['f'] = local_f
    polynomials['g'] = local_g
    f.set(polynomials['f'].Get_polystr())
    g.set(polynomials['g'].Get_polystr())


def set_f():
    polynomials['f'] = Poly(String_to_Poly(f.get()))
    polynomials['f'].Get_polystr()


def set_g():
    polynomials['g'] = Poly(String_to_Poly(g.get()))
    polynomials['g'].Get_polystr()


def Multiply():
    f = polynomials['f']

    g = polynomials['g']

    h = f.Multiply(g)

    polynomials['h'] = h
    results.set(h.Get_polystr())


def Add():
    f = polynomials['f']
    f.Get_polystr()
    g = polynomials['g']
    g.Get_polystr()

    h = f.Add(g)

    results.set(h.Get_polystr())

def Simplify():
    print('Simplifying~~~~~~~~~~~~~~')
    local_f = polynomials['f']
    local_g = polynomials['g']
    local_h = Poly(String_to_Poly(results.get()))

    local_f.Simplify()
    print('simplified f')
    local_g.Simplify()
    print('simplified g')
    local_h.Simplify()
    print('simplified h')

    results.set(local_h.Get_polystr())
    f.set(local_f.Get_polystr())
    g.set(local_g.Get_polystr())

    polynomials['f'] = local_f
    polynomials['g'] = local_g

def Graph_Points():

    X = dictionaries['values']

    h = Poly(String_to_Poly(results.get()))

    h.Get_polystr()
    Y = h.Eval_Points(points=X)
    Y_labels = []
    for y in Y:
        Y_labels.append(dictionaries['dictionary'][str(y)])
    print(Y)
    print(Y_labels)

    fig, ax1 = plt.subplots()

    ax2 = ax1.twinx()
    ax1.plot(X, Y, 'g-')

    ax2.plot(X, Y_labels, 'b-')

    ax1.set_xlabel('X data')
    ax1.set_ylabel('Y1 data', color='g')
    ax2.set_ylabel('Y2 data', color='b')

    #plt.show()
    plt.savefig('graph.png')
    plt.clf()

    # plt.scatter(X,Y)
    # plt.savefig('graph.png')
    # plt.clf()



    img = (Image.open("graph.png"))
    image = ImageTk.PhotoImage(img.resize((800, 600), Image.ANTIALIAS))
    imageLabel.configure(image=image)
    imageLabel.image = image



if __name__ == "__main__":
    sand_window = tk.Tk()
    sand_window.attributes('-fullscreen', True)
    container = tk.Frame(sand_window)
    canvas = tk.Canvas(container)
    scrollbarh = tk.Scrollbar(container, orient="horizontal", command=canvas.xview)
    scrollbarv = tk.Scrollbar(container, orient="vertical", command=canvas.yview)
    canvas.configure(yscrollcommand=scrollbarv.set)
    canvas.configure(xscrollcommand=scrollbarh.set)
    scrollbarh.pack(side="bottom", fill="y")
    scrollbarv.pack(side="right", fill="y")
    container.pack(fill=tk.BOTH, expand=True)
    canvas.pack(side="left", fill="both", expand=True)
    scrollable_frame = tk.Frame(canvas)
    scrollable_frame.bind(
        "<Configure>",
        lambda e: canvas.configure(
            scrollregion=canvas.bbox("all")
        )
    )
    canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")

    label_poly = tk.StringVar()
    button_quit = tk.Button(scrollable_frame, text="QUIT",
                            command=close)
    button_quit.grid(row=0, column=0)

    label_poly.set("Type polynomials in increasing order of coefficient powers and variable powers")
    label_poly_Dir = tk.Label(scrollable_frame, textvariable=label_poly)
    label_poly_Dir.grid(row=0, column=1)

    label_field_size = tk.StringVar()
    label_field_size.set("Field Size:")
    label_field_size_Dir = tk.Label(scrollable_frame, textvariable=label_field_size)
    label_field_size_Dir.grid(row=0, column=3)
    field_size_tk = tk.IntVar()
    field_size_tk.set(parameters['field_size'])
    field_size_entry = tk.Entry(scrollable_frame, textvariable=field_size_tk)
    field_size_entry.insert(tk.END, '')
    field_size_entry.grid(row=0, column=4)

    # base_tk.trace('w', update_functions)
    button_update_functions = tk.Button(scrollable_frame, text="Update Functions",
                                        command=update_functions)
    button_update_functions.grid(row=0, column=5)
    label_f = tk.StringVar()
    label_f.set("Polynomial f:")
    label_f_Dir = tk.Label(scrollable_frame, textvariable=label_f)
    label_f_Dir.grid(row=1, column=0)
    f = tk.StringVar()
    f_entry = tk.Entry(scrollable_frame, textvariable=f)
    f_entry.insert(tk.END, '')
    f_entry.grid(row=2, column=0)
    #f.trace('w', update_function_f)
    button_set_f = tk.Button(scrollable_frame, text="Set function f",
                               command=set_f)
    button_set_f.grid(row=3, column=0)
    label_g = tk.StringVar()
    label_g.set("Polynomial g:")
    label_g_Dir = tk.Label(scrollable_frame, textvariable=label_g)
    label_g_Dir.grid(row=1, column=1)
    g = tk.StringVar()
    g_entry = tk.Entry(scrollable_frame, textvariable=g)
    g_entry.insert(tk.END, '')
    g_entry.grid(row=2, column=1)
    #g.trace('w', update_function_g)
    button_randomize = tk.Button(scrollable_frame, text="Randomize",
                                 command=Randomize)
    button_randomize.grid(row=2, column=5)
    button_set_g = tk.Button(scrollable_frame, text="Set Function g",
                               command=set_g)
    button_set_g.grid(row=3, column=1)

    button_graph = tk.Button(scrollable_frame,
                             text="GRAPH results",command=Graph_Points)
    button_graph.grid(row=4, column=0)

    button_multiply = tk.Button(scrollable_frame, text="Multiply",
                                command=Multiply)
    button_multiply.grid(row=4, column=1)
    button_add = tk.Button(scrollable_frame, text="Add",
                           command=Add)
    button_add.grid(row=4, column=2)

    button_simplify = tk.Button(scrollable_frame,text="Simplify",command=Simplify)
    button_simplify.grid(row=4,column=3)
    label_results = tk.StringVar()
    label_results.set("Results")
    label_results_Dir = tk.Label(scrollable_frame, textvariable=label_results)
    label_results_Dir.grid(row=5, column=0)
    # g = tk.StringVar()
    # g_entry = tk.Entry(scrollable_frame, textvariable=g)
    # g_entry.insert(tk.END, '')
    # g_entry.grid(row=2, column=1)
    results = tk.StringVar()
    results_value = tk.Entry(scrollable_frame,textvariable = results)#tk.Label(scrollable_frame, textvariable=results)
    results_value.insert(tk.END,'')
    results_value.grid(row=5, column=1)
    imageLabel = tk.Label(scrollable_frame)
    imageLabel.grid(row=6, column=0,columnspan=4)
    sand_window.mainloop()
