class WildParticle:
    def __init__(self, x, y, u, v, mass, color, lifetime):
        self.fps = 24
        self.position = { 'x': x, 'y': y }
        self.velocity = { 'u': float(u) / self.fps, 'v': float(v) / self.fps }
        self.acceleration = { 'a': 0, 'b': 0 }
        self.mass = mass
        self.color = color
        self.death = self.fps * lifetime
        self.time = 0
        print("Added particle at ({}, {}): mass = {}; speed: ({}, {})".format(
            self.position['x'], self.position['y'], self.mass,
            self.velocity['u'], self.velocity['v']))

    def __str__(self):
        return "Particle at ({}, {}): m = {}; speed: ({}, {})".format(
            self.position['x'], self.position['y'], self.mass,
            self.velocity['u'], self.velocity['v'])

    __repr__ = __str__
