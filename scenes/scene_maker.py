import random

num_particles = 3000
bbox_mins = [-2, -2, -2]
bbox_maxs = [2, 2, 2]
pad = [0.5, 0.5, 0.5]
particles = []

for i in range(num_particles):
    particles.append((
        random.uniform(bbox_mins[0]+pad[0], bbox_maxs[0]-pad[0]),
        random.uniform(bbox_mins[1]+pad[1], bbox_maxs[1]-pad[1]),
        random.uniform(bbox_mins[2]+pad[2], bbox_maxs[2]-pad[2]),
    ))

f = open("3000_random.txt", "w")
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
