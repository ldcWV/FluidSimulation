import random

num_particles = 10000
bbox_mins = [-0.5, 0, -0.5]
bbox_maxs = [0.5, 8, 0.5]
x_interval = (-0.5, 0)
y_interval = (1, 4)
z_interval = (-0.25, 0.25)
particles = []

for i in range(num_particles):
    particles.append((
        random.uniform(x_interval[0], x_interval[1]),
        random.uniform(y_interval[0], y_interval[1]),
        random.uniform(z_interval[0], z_interval[1])
    ))

f = open(f"{num_particles}_small.txt", "w")
f.write(str(num_particles) + "\n")
for i in range(3):
    f.write(str(bbox_mins[i]) + " ")
for i in range(3):
    f.write(str(bbox_maxs[i]))
    if i != 2:
        f.write(" ")
f.write("\n")
for i, (x,y,z) in enumerate(particles):
    f.write(str(i) + " " + str(x) + " " + str(y) + " " + str(z) + " 0.0 0.0 0.0 0.0 0.0 0.0\n")
f.close()
