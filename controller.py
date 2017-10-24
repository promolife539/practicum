import math
import model
import cython_verlet
import threading
import numpy
import matplotlib
import matplotlib.pyplot as plot
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg


def pure_python_verlet(particles):
    G = 6.67 * (10 ** -11)
    for p in particles:
        new_a = sum(map(lambda e:
            G * e.mass * (e.position['x'] - p.position['x'])
            / (math.sqrt(
                (e.position['x'] - p.position['x'])**2 +
                (e.position['y'] - p.position['y'])**2
            )) , filter(lambda e: not e is p, particles)))
        new_b = sum(map(lambda e:
            G * e.mass * (e.position['y'] - p.position['y'])
            / (math.sqrt(
                (e.position['x'] - p.position['x'])**2 +
                (e.position['y'] - p.position['y'])**2
            )), filter(lambda e: not e is p, particles)))
        if p.time > 0:
            p.position['x'] += p.velocity['u'] + 0.5 * p.acceleration['a']
            p.position['y'] += p.velocity['v'] + 0.5 * p.acceleration['b']
        p.time += 1
        p.velocity['u'] += 0.5 * (new_a + p.acceleration['a'])
        p.velocity['v'] += 0.5 * (new_b + p.acceleration['b'])
        p.acceleration['a'] = new_a
        p.acceleration['b'] = new_b
    return [p.position for p in particles]

def verlet_worker(particles, begin, end):
    G = 6.67 * (10 ** -11)
    for p in particles[begin : end]:
        new_a = sum(map(lambda e:
            G * e.mass * (e.position['x'] - p.position['x'])
            / (math.sqrt(
                (e.position['x'] - p.position['x'])**2 +
                (e.position['y'] - p.position['y'])**2
            )) , filter(lambda e: not e is p, particles)))
        new_b = sum(map(lambda e:
            G * e.mass * (e.position['y'] - p.position['y'])
            / (math.sqrt(
                (e.position['x'] - p.position['x'])**2 +
                (e.position['y'] - p.position['y'])**2
            )), filter(lambda e: not e is p, particles)))
        if p.time > 0:
            p.position['x'] += p.velocity['u'] + 0.5 * p.acceleration['a']
            p.position['y'] += p.velocity['v'] + 0.5 * p.acceleration['b']
        p.time += 1
        p.velocity['u'] += 0.5 * (new_a + p.acceleration['a'])
        p.velocity['v'] += 0.5 * (new_b + p.acceleration['b'])
        p.acceleration['a'] = new_a
        p.acceleration['b'] = new_b

def multiprocess_verlet(particles):
    jobs = []
    for i in range(len(particles)):
        job = threading.Thread(target = verlet_worker, args = (particles, i, i + 1))
        job.start()
        jobs.append(job)
    for j in jobs:
        j.join()
    return [p.position for p in particles]

class ParticlePlot(FigureCanvasQTAgg):
    def __init__(self, parent, width, height, dpi, size_policy):
        figure = matplotlib.figure.Figure(figsize = (width, height), dpi = dpi,
            facecolor = 'white')
        self.axes = figure.add_axes([0.005,0.005,0.990,0.990], frameon=True, aspect=1)
        FigureCanvasQTAgg.__init__(self, figure)
        self.setParent(parent)
        FigureCanvasQTAgg.setSizePolicy(self, size_policy, size_policy)
        FigureCanvasQTAgg.updateGeometry(self)
        self.figure.canvas.draw()

    def update_plot(self, particles, updater):
        self.axes.cla()
        self.axes.set_xlim(-45, 45), self.axes.set_xticks([])
        self.axes.set_ylim(-25, 25), self.axes.set_yticks([])
        data = updater(particles)
        mass_data = [ p.mass / 1000 for p in particles ]
        color_data = [ p.color for p in particles ]
        x_data = [ p['x'] for p in data ]
        y_data = [ p['y'] for p in data ]
        self.scatter = self.axes.scatter(x_data, y_data, s = mass_data, lw = 0.5,
            c = color_data)
        self.figure.canvas.draw()

class ParticleController:
    defaults = {
        'mass': 1250 * 1000,
        'lifetime': 4,
        'velocity': { 'u': 5, 'v': 7 },
        'position': { 'x': 0, 'y': 0 },
        'color': (0, 1, 0),
        'method': 0
    }

    def __init__(self):
        self.__mass = __class__.defaults['mass']
        self.__lifetime = __class__.defaults['lifetime']
        self.__velocity = __class__.defaults['velocity']
        self.__position = __class__.defaults['position']
        self.__color = __class__.defaults['color'];
        self.method = __class__.defaults['method']
        self.particles = []
        self.updaters = [
            pure_python_verlet,
            cython_verlet.cython_verlet,
            multiprocess_verlet,
        ]
        self.methods = [
            "Pure Python Verlet algorithm implementation",
            "Cython Verlet algorithm implementation",
            "Multiprocess Verlet algorithm implementation",
        ]

    def __add_particle(self):
        self.particles.append(model.WildParticle(
            self.position['x'], self.position['y'],
            self.velocity['u'], self.velocity['v'],
            self.mass, self.color, self.lifetime
        ) )

    @property
    def position(self):
        return self.__position

    @position.setter
    def position(self, value):
        self.__position = value

    @property
    def velocity(self):
        return self.__velocity

    @velocity.setter
    def velocity(self, value):
        self.__velocity = value

    @property
    def acceleration(self):
        return self.__acceleration

    @acceleration.setter
    def acceleration(self, value):
        self.__acceleration = value

    @property
    def mass(self):
        return self.__mass

    @mass.setter
    def mass(self, value):
        self.__mass = value

    @property
    def lifetime(self):
        return self.__lifetime

    @lifetime.setter
    def lifetime(self, value):
        self.__lifetime = value

    @property
    def color(self):
        return self.__color

    @color.setter
    def color_set(self, value):
        self.__color = value
