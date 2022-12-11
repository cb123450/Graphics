import tkinter as tk
from tkinter import ttk
from tkinter.messagebox import showinfo
from tkinter import Canvas
from tkinter import *

#create window
window = tk.Tk()
frame = Frame(window)
frame.pack()
window.title("Mini-project")


#window dimensions
window.configure(width=1500, height=800)
window.configure(bg='lightgray')

# move window to the center of the screen
windowWidth = window.winfo_reqwidth()
windowHeight = window.winfo_reqheight()
posRight = int(window.winfo_screenwidth() / 2 - windowWidth / 2)
posDown = int(window.winfo_screenheight() / 2 - windowHeight / 2)
window.geometry("+{}+{}".format(posRight, posDown))

#create canvas to draw and bind it to the window
canvas = Canvas(window, height=windowHeight-100, width=windowWidth-100)
canvas.pack()

global x, y
global radius
radius = 10

vertex_centers = []
circles = []
edge_vertices = []
lines = []

def place_vertex(event):

    if (vertex_bool.get() == 1 and edge_bool.get() == 0 and color_bool.get() == 0 and blue_bool.get() == 0):
        # x increases from left to right
        # y increases from up to down
        x0 = event.x - radius
        y0 = event.y - radius
        x1 = event.x + radius
        y1 = event.y + radius
        center = ((x0+x1)/2.0, (y0+y1)/2.0)
        vertex_centers.append(center)

        tag = (center[0], center[1])

        #fill the circle black
        circle = canvas.create_oval(x0, y0, x1, y1, fill = "black", tags = tag)
        circles.append(circle)

        return 0

#helper function to determine if the user clicked on a vertex already drawn on the screen
def click_on_vert(click_x, click_y):
    for center in vertex_centers:
        #bounding box for a given vertex
        left = center[0]-radius
        right = center[0]+radius
        top = center[1]-radius
        bottom = center[1]+radius

        if (left <= click_x and click_x <= right and top <= click_y and click_y <= bottom):
            return center

    #if vertex not found, no vertex was clicked on so return 0         
    return 0
    

def draw_edge(event):  
    if (vertex_bool.get() == 0 and edge_bool.get() == 1 and color_bool.get() == 0 and blue_bool.get() == 0):
        ret = click_on_vert(event.x, event.y)
        if (ret != 0):
            edge_vertices.append(ret)

        if (len(edge_vertices) > 0 and len(edge_vertices)%2 ==0):
            last = edge_vertices[len(edge_vertices)-1]
            second_to_last = edge_vertices[len(edge_vertices)-2]
            tag = (last[0], last[1], second_to_last[0], second_to_last[1])
            line = canvas.create_line(last[0], last[1], second_to_last[0], second_to_last[1], tags = tag)
            lines.append(line)

#CHANGE THE COLORS OF CIRCLES
def color_vertex_blue(event):
    if (vertex_bool.get() == 0 and edge_bool.get() == 0 and color_bool.get() ==1 and blue_bool.get() == 1 and red_bool.get() == 0 and green_bool.get() == 0 and yellow_bool.get() == 0 and black_bool.get() == 0 and purple_bool.get() == 0 and pink_bool.get() == 0):
        ret = click_on_vert(event.x, event.y)
        if (ret != 0):
            for circle in circles:
                tag = canvas.gettags(circle)
                if ((abs(float(ret[0]) - float(tag[0]))) < .0001 and (abs(float(ret[1]) - float(tag[1]))) < .0001):
                    canvas.itemconfig(circle, fill = "blue")

def color_vertex_red(event):
    if (vertex_bool.get() == 0 and edge_bool.get() == 0 and color_bool.get() ==1 and blue_bool.get() == 0 and red_bool.get() == 1 and green_bool.get() == 0 and yellow_bool.get() == 0 and black_bool.get() == 0 and purple_bool.get() == 0 and pink_bool.get() == 0):
        ret = click_on_vert(event.x, event.y)
        if (ret != 0):
            for circle in circles:
                tag = canvas.gettags(circle)
                if ((abs(float(ret[0]) - float(tag[0]))) < .0001 and (abs(float(ret[1]) - float(tag[1]))) < .0001):
                    canvas.itemconfig(circle, fill = "red")


#variables
vertex_bool = tk.IntVar()
edge_bool = tk.IntVar()
color_bool = tk.IntVar()
blue_bool = tk.IntVar()
red_bool = tk.IntVar()
green_bool = tk.IntVar()
yellow_bool = tk.IntVar()
black_bool = tk.IntVar()
purple_bool = tk.IntVar()
pink_bool = tk.IntVar()

def place_vertex_or_edge():
    if (vertex_bool.get() == 1 and edge_bool.get() == 0 and color_bool.get() == 0 and blue_bool.get() == 0):
        #vertex mode
        window.bind('<Control-Button-1>', place_vertex)
        
    if (vertex_bool.get() == 0 and edge_bool.get() == 1 and color_bool.get() == 0 and blue_bool.get() == 0):
        #edge mode
        window.bind('<Control-Button-1>', draw_edge)

def color_mode():
    if (vertex_bool.get() == 0 and edge_bool.get() == 0 and color_bool.get() ==1): 
        if (blue_bool.get() == 1 and red_bool.get() == 0 and green_bool.get() == 0 and yellow_bool.get() == 0 and black_bool.get() == 0 and purple_bool.get() == 0 and pink_bool.get() == 0):
            #color blue
            window.bind('<Button-1>', color_vertex_blue)
        if (blue_bool.get() == 0 and red_bool.get() == 1 and green_bool.get() == 0 and yellow_bool.get() == 0 and black_bool.get() == 0 and purple_bool.get() == 0 and pink_bool.get() == 0):
            #color red
            window.bind('<Button-1>', color_vertex_red)
        

vertexCheckButton = Checkbutton(frame, text='Vertex mode: control click to place vertices', command=place_vertex_or_edge, variable=vertex_bool)
vertexCheckButton.pack(side = BOTTOM)

edgeCheckButton = Checkbutton(frame, text='Edge mode: control click to place edge', command=place_vertex_or_edge, variable=edge_bool)
edgeCheckButton.pack(side = BOTTOM)

colorCheckButton = Checkbutton(frame, text='Color mode: turn on to color vertices', command=color_mode, variable=color_bool)
colorCheckButton.pack(side = BOTTOM)

blueCheckbutton = Checkbutton(frame, text='Blue: left click to color vertices blue', command=color_mode, variable=blue_bool, fg = "blue")
blueCheckbutton.pack(side = BOTTOM)

redCheckbutton = Checkbutton(frame, text='Red: left click to color vertices red', command=color_mode, variable=red_bool, fg = "red")
redCheckbutton.pack(side = BOTTOM)

greenCheckbutton = Checkbutton(frame, text='Green: left click to color vertices green', command=color_mode, variable=green_bool, fg = "green")
greenCheckbutton.pack(side = BOTTOM)

yellowCheckbutton = Checkbutton(frame, text='Yellow: left click to color vertices yellow', command=color_mode, variable=yellow_bool, fg = "yellow")
yellowCheckbutton.pack(side = BOTTOM)

blackCheckbutton = Checkbutton(frame, text='Black: left click to color vertices black', command=color_mode, variable=black_bool, fg = "black")
blackCheckbutton.pack(side = BOTTOM)

purpleCheckbutton = Checkbutton(frame, text='Purple: left click to color vertices purple', command=color_mode, variable=purple_bool, fg = "purple")
purpleCheckbutton.pack(side = BOTTOM)

pinkCheckbutton = Checkbutton(frame, text='Pink: left click to color vertices pink', command=color_mode, variable=pink_bool, fg = "pink")
pinkCheckbutton.pack(side = BOTTOM)

window.bind(place_vertex_or_edge)

window.mainloop()