
import numpy as np
import matplotlib.pyplot as plt

# run simulation
p_mit = 80  # probability of dividing
p_death = 25  # probability of dying
time_clear = 3  # time to clear after death
p_symm = 50  # probability that division is symmetrical
d_max = 15

T = 15 # number of iteraitons
N = 10  # size of grid
Q = 2  # number of cycles
time = np.zeros((N, N)) # grid to keep track of time before clearing grid

#function returns new array after applying probabilities of cell death/division etc
def advance(grid, N, time_clear, p_mit, p_death, time):
    result = grid.copy()
    count3 = 0
    count2 = 0
    count1 = 0
    count0 = 0
    # loop over all cells in grid (except border cells)
    for i in range(1, N - 1):
        for j in range(1, N - 1):
            neighbours = np.array([[i + 1, j], [i - 1, j], [i, j + 1], [i, j - 1]])
            # if cell i,j is a progenitor cell and at least one neighbouring cell is empty
            if grid[i, j] == 1:
                count1 += 1
                if np.random.choice(100) <= p_mit:
                    # symmetric division
                    if result[i + 1, j] == 0 or result[i - 1, j] == 0 or result[i, j + 1] == 0 or result[i, j - 1] == 0:
                        # randomly select an empty neighbouring cell and set it 1 (i.e. progenitor)
                        # NB there may be more efficient ways to do this
                        hit = False
                        while hit == False:
                            r = np.random.choice(4)
                            x, y = neighbours[r, :]
                            if result[x, y] == 0:
                                result[x, y] = 1
                                hit = True
                        if np.random.choice(100) <= p_death:
                            result[x, y] = 2

                            # d_max
                            count_prog = 0
                            for k in range(N):
                                for m in range(N):
                                    if result[k, m] == 1:
                                        count_prog += 1
                            if count_prog >= d_max:
                                result[i, j] = 4
                                result[x, y] = 2
                        else:
                            count_prog = 0
                            for k in range(N):
                                for m in range(N):
                                    if result[k, m] == 1:
                                        count_prog += 1

                            if count_prog >= d_max:
                                result[x, y] = 4

            # if stem cell
            if grid[i, j] == 3:
                count3 += 1
                if np.random.choice(100) <= p_mit:

                    # symmetric division
                    if result[i + 1, j] == 0 or result[i - 1, j] == 0 or result[i, j + 1] == 0 or result[i, j - 1] == 0:
                        if np.random.choice(100) <= p_symm:

                            hit = False
                            while hit == False:
                                r = np.random.choice(4)
                                x, y = neighbours[r, :]
                                if result[x, y] == 0:
                                    result[x, y] = 3
                                    hit = True

                        # asymmetric division
                        else:

                            # asymmetric, 1 dead
                            if np.random.choice(100) <= p_death:
                                hit = False
                                while hit == False:
                                    r = np.random.choice(4)
                                    x, y = neighbours[r, :]
                                    if result[x, y] == 0:
                                        result[x, y] = 2
                                        hit = True

                            else:
                                # asymmetric, progenitor
                                hit = False
                                while hit == False:
                                    r = np.random.choice(4)
                                    x, y = neighbours[r, :]
                                    if result[x, y] == 0:
                                        result[x, y] = 1
                                        hit = True

    # loops over all cells (because even border cells can die)
    for i in range(N):
        for j in range(N):
            # if alive, maybe die
            # kills off too many cells!
            # if grid[i,j] == 1 or grid[i,j] == 3:
            # if np.random.choice(100)<=p_death:
            # result[i,j] = 2 # if r = 0, kill the cell

            # if dead, count how long its been dead for
            if grid[i, j] == 2:
                count2 += 1
                # 2 means  dead, 1 means alive , 0 means empty:limbo?
                time[i, j] += 1
                # if dead too long, empty
                if time[i, j] == time_clear:
                    time[i, j] = 0
                    result[i, j] = 0
                else:
                    result[i, j] = 2

    return result, time, count1, count2, count3,


# returns it as a colourful grid rather than an array of numbers
automaton = np.zeros((T, N, N))
x = np.zeros((N, N))

# x[np.random.choice(N),np.random.choice(N)] = 1 # sets random coordinate as prog cell can change this to try different starting positions,
# x[np.random.choice(N),np.random.choice(N)] = 3 # sets random coordinate as stem cell

# starting each in the middle
x[int(N / 2) + 1, int(N / 2) + 1] = 1  # progenitor celliin middle
x[int(N/2) - 1,int(N/2) - 1] = 3 # stem cell in middle

# x[int(N/2)+1,int(N/2)+1] = 1 # sets progenitor cell starting point
# x[int(N/2)-1,int(N/2)-1] = 3 # sets stem cell starting point

count = np.zeros((T, 3))
plt.figure(figsize=(12, 2))
plt.title("Visual of Neurospheres")
plt.axis('off')
for i in range(T):
    # print(x)
    automaton[i, :, :] = x
    x, t, a, b, c, = advance(x, N, time_clear, p_mit, p_death, time)
    count[i, :] = a, b, c
    # print(a,c)
    plt.subplot(1, T, i + 1)
    plt.imshow(automaton[i, :, :], cmap="Accent", vmin=0, vmax=4)
    plt.axis('off')
    plt.savefig('Visualisation of neurospheres.pdf')


# plot number of each cells over time (1 unit of time = 1 iteration)
plt.figure(figsize=(8, 5))
plt.plot(count[:, 0], label='Progenitor cells', color='orange')  # first column
plt.plot(count[:, 1], label='Dead cells', color='blue')  # second column
plt.plot(count[:, 2], label='Stem cells', color='red')
plt.xlabel("Iterations")  # add an x-axis label
plt.ylabel("Cell Count #")  # add a y-axis label
plt.legend()
plt.title("Cell numbers over time")
plt.savefig('Cell numbers over time.pdf')
plt.show()

# meaning of colours
# green: empty, 0
# light orange: progenitor, 1
# blue: dead, 2
# brown: stem cell, 3
# grey: differentiated, 4