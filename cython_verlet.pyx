import math

def cython_verlet(particles):
    G_num = 0.0667
    G_den = 1000000000
    for p in particles:
        new_a = sum(map(lambda e:
            G_num * e.mass * (e.position['x'] - p.position['x'])
            / (math.sqrt(
                (e.position['x'] - p.position['x'])**2 +
                (e.position['y'] - p.position['y'])**2
            ) * G_den) , filter(lambda e: not e is p, particles)))
        new_b = sum(map(lambda e:
            G_num * e.mass * (e.position['y'] - p.position['y'])
            / (math.sqrt(
                (e.position['x'] - p.position['x'])**2 +
                (e.position['y'] - p.position['y'])**2
            ) * G_den), filter(lambda e: not e is p, particles)))
        if p.time > 0:
            p.position['x'] += p.velocity['u'] + 0.5 * p.acceleration['a']
            p.position['y'] += p.velocity['v'] + 0.5 * p.acceleration['b']
        p.time += 1
        p.velocity['u'] += 0.5 * (new_a + p.acceleration['a'])
        p.velocity['v'] += 0.5 * (new_b + p.acceleration['b'])
        p.acceleration['a'] = new_a
        p.acceleration['b'] = new_b
    return [p.position for p in particles]
