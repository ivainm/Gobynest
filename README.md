# Goby Nest data analysis

## Script files:

### nestcode.R 

Code to produce all the statistical analysis and figures of the mansucript, except the path analysis figure (handmade).

## Data files

### Coord.txt

Coordinates of the different field populations.\
"Area": name of the field population\
"Lat" and "Long": Latitude and Longitude of the field population

### eggsize.txt

Data of egg diameter.\
"Area": name of the field population\
"Location": name of field site within the population\
"Timeseason": period of the breeding season, either peak or late\
"NestID": Nest identifier, unique within Location\
"Nest_size": either large or small\
"dx.y": diameter of egg x, replicate y in pixels\
"scale": pixel measurement of scale bar\
"eyed": egg developmental stage, were eyes visible?

### env_info.txt

Thermal environment data for each field population.\
"Area": name of the field population\
"Broods": estimated number of broods that could be produced consecutively over the breeding period\
"Ddays": number of degree days accumulated over the breeding period\
"Ddaysfixed": number of degree days over a fixed time window\
"Broodsfixed": number of consecutive broods over a fixed time window

### fullnest2.txt

Data collected from the nests.\
"area": name of the field population\
"location": name of field site within the population\
"Timepoint": period of the breeding season, either peak (2) or late (3)\
"date": collection date\
"NestID": Nest identifier, unique within Location\
"Nest_size": either large or small\
"Live_eggs_field": Field notes, are there life eggs in the nest?\
"Egg_squares_field": Field notes, estimates the surface covered by eggs, in grid squares\
"Eyed_field": Field notes, are there eyed eggs?\
"male": Field notes, was a guarding male captured in the nest?\
"length" "weight": male measurements in mm and centigrams\
"comment_analysis": comment made during the picture analysis stage\
"Eyed_analysis": Picture analysis, are there eyed eggs?\
"Live_eggs_analysis": Picture analysis, are there live eggs?\
"TotBrood_surface": Total brood surface in pixels\
"1sq_surface": Surface in pixel of 1 square of the grid pattern\
"Cann_surface": Cannibalised surface in pixels\
"Grown_surface": Overgrown surface in pixels\
"Egg_surface": Surface in pixels of the square mask used to count eggs\
"Egg1",	"Egg2",	"Egg3",	"Egg4",	"Egg5": 5 replicated egg counts using the square mask\
"square_length": Number of grid squares in the length of the nest inside surface\
"square_height": Number of grid squares in the width/height of the nest inside surface

### testpanel.txt

Contains the identity of the nests used for nest analysis method development (supplementary file F2 of the manuscript)
