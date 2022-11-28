import matplotlib.pyplot as plt

from copy import deepcopy
import random
import math


def generate_cubic_lattice(num_atoms, density):
    """
    Generate points on a cubic lattice using a desired final density.

    Parameters
    ----------
    num_atoms: int
      The number of atoms to place on the lattice.
    density: float
      The desired system density.

    Returns
    -------
    coords: list
      A nested list of generated coordinates.
    """

    # Calculate box length based on number of atoms and density.
    volume = num_atoms / density
    box_length = math.pow(volume, (1.0 / 3.0))

    # Calculate the upper bound of cube size.
    # Our approach will be to place atoms until
    # we place all needed. For this, we need
    # to determine the maximum number of atoms on each
    # side.
    max_side = math.ceil(math.pow(num_atoms, (1.0 / 3.0)))

    # Determine spacing based on number of atoms
    # and box length.
    spacing = box_length / max_side  # units length / atom

    coordinates = []
    count = 0

    for i in range(max_side):
        for j in range(max_side):
            for k in range(max_side):
                coordinates.append([i * spacing, j * spacing, k * spacing])
                count += 1
                if count == num_atoms:
                    return coordinates, box_length


def calculate_distance_vector(coord1, coord2, box_length=None):
    """
    Calculate the distance vector between two points.
    """

    dist_vec = []

    for i in range(3):
        dim_dist = coord1[i] - coord2[i]

        if box_length:
            dim_dist = dim_dist - box_length * round(dim_dist / box_length)

        dist_vec.append(dim_dist)

    return dist_vec


def calculate_distance(coord1, coord2, box_length=None):
    """
    Calculate the distance between two 3D coordinates.

    Parameters
    ----------
    coord1, coord2: list
        The atomic coordinates

    Returns
    -------
    distance: float
        The distance between the two points.
    """

    dist_vect = calculate_distance_vector(coord1, coord2, box_length)
    distance = math.sqrt(sum([x**2 for x in dist_vect]))
    return distance


def calculate_LJ(r_ij):
    """
    The LJ interaction energy between two particles.

    Computes the pairwise Lennard Jones interaction energy based
    on the separation distance in reduced units.

    Parameters
    ----------
    r_ij : float
        The distance between the particles in reduced units.

    Returns
    -------
    pairwise_energy : float
        The pairwise Lennard Jones interaction energy in reduced units.

    """

    r6_term = math.pow(1 / r_ij, 6)
    r12_term = math.pow(r6_term, 2)

    pairwise_energy = 4 * (r12_term - r6_term)

    return pairwise_energy


def calculate_force_vector(distance_vector):
    """Calculate the force vector based on a distance vector and the LJ potential"""

    r2 = sum([i**2 for i in distance_vector])
    f = [0, 0, 0]
    inv = 1 / r2

    for d in range(3):
        r6_term = inv**3
        ff = 48 * inv * r6_term * (r6_term - 0.5)
        f[d] = ff * distance_vector[d]
    return f


def calculate_forces(coordinates, box_length, cutoff):
    """Calculate forces on all particles"""

    n_particles = len(coordinates)

    forces = []
    for i in range(n_particles):
        forces.append([0, 0, 0])

    for i in range(n_particles):
        for j in range(i + 1, n_particles):

            dist_vector = calculate_distance_vector(
                coordinates[i], coordinates[j], box_length
            )

            r2 = sum([i**2 for i in dist_vector])

            if r2 < cutoff**2:
                force_vector = calculate_force_vector(dist_vector)

                for k in range(len(force_vector)):
                    forces[i][k] += force_vector[k]
                    forces[j][k] -= force_vector[k]

    return forces


def half_step(positions, velocities, forces, timestep):
    """Perform half a step of a Velocity Verlet integration"""

    positions = deepcopy(positions)

    n_particles = len(positions)

    half_velocities = []

    for i in range(n_particles):
        # calculate half velocities
        v_half_x = velocities[i][0] + 0.5 * forces[i][0] * timestep
        v_half_y = velocities[i][1] + 0.5 * forces[i][1] * timestep
        v_half_z = velocities[i][2] + 0.5 * forces[i][2] * timestep

        # calculate new positions
        positions[i][0] += v_half_x * timestep
        positions[i][1] += v_half_y * timestep
        positions[i][2] += v_half_z * timestep

        half_velocities.append([v_half_x, v_half_y, v_half_z])

    return positions, half_velocities


def set_initial_velocities(temperature, num_atoms):
    """Set initial velocities based on a target temperature and number of atoms."""

    velocities = []
    velocity_sum = [0, 0, 0]
    velocity_sum_sq = 0

    for i in range(num_atoms):
        atom_velocity = []
        for j in range(3):
            v_j = random.random() - 0.5
            velocity_sum[j] += v_j
            velocity_sum_sq += v_j**2

            atom_velocity.append(v_j)

        velocities.append(atom_velocity)

    # velocity center of mass
    velocity_sum = [x / num_atoms for x in velocity_sum]

    # mean-squared velocity
    velocity_sum_sq = velocity_sum_sq / num_atoms

    scale_factor = math.sqrt(3 * temperature / velocity_sum_sq)

    for i in range(num_atoms):
        for j in range(3):
            velocities[i][j] = (velocities[i][j] - velocity_sum[j]) * scale_factor

    return velocities


def calculate_temperature(velocities):
    """Calculate temperature based on velocities."""

    n_particles = len(velocities)

    temperature = 0
    for i in range(n_particles):
        v_s = sum([x**2 for x in velocities[i]])
        temperature += v_s / (3 * n_particles - 1)

    return temperature


def calculate_potential_energy(coordinates, box_length, cutoff):
    """
    Calculate the total Lennard Jones energy of a system of particles.

    Parameters
    ----------
    coordinates : list
        Nested list containing particle coordinates.

    Returns
    -------
    total_energy : float
        The total pairwise Lennard Jones energy of the system of particles.
    """

    total_energy = 0

    num_atoms = len(coordinates)

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):

            dist_ij = calculate_distance(
                coordinates[i], coordinates[j], box_length=box_length
            )

            if dist_ij < cutoff:
                interaction_energy = calculate_LJ(dist_ij)
                total_energy += interaction_energy

    return total_energy


def write_xyz(coordinates, fn, mode="w+"):

    with open(fn, mode) as f:
        f.write(f"{len(coordinates)}\n")
        f.write("#comment\n")

        for i in range(len(coordinates)):
            f.write(f"H {coordinates[i][0]} {coordinates[i][1]} {coordinates[i][2]}\n")


if __name__ == "__main__":

    n_particles = 125
    target_temp = 0.9
    target_density = 0.9
    num_steps = 1000
    cut = 3
    cut2 = cut**2
    timestep = 0.01
    total_energies = []

    velocities = set_initial_velocities(target_temp, n_particles)

    positions, box_length = generate_cubic_lattice(n_particles, target_density)

    initial_pe = calculate_potential_energy(positions, box_length, cut)

    temperature = calculate_temperature(velocities)

    print(f"Target temperature set to {target_temp}.")
    print(f"Temperature calculated to {temperature}")
    print(f"Initial potential energy {initial_pe/n_particles}")
    print(f"The box length is {box_length}")
    print(f"The density is {n_particles/box_length**3}")

    write_xyz(positions, f"trajectory_{n_particles}_{target_density}.xyz")

    for s in range(num_steps):

        # MD simulation loop.

        # calculate forces on particles
        forces = calculate_forces(positions, box_length, cut)

        # integrate equations of motion (apply half of Velocity Verlet)
        positions, half_velocities = half_step(positions, velocities, forces, timestep)

        # other half of velocity verlet
        forces = calculate_forces(
            positions, box_length, cut
        )  # need forces with updated positions

        for i in range(n_particles):
            velocities[i][0] = half_velocities[i][0] + 0.5 * forces[i][0] * timestep
            velocities[i][1] = half_velocities[i][1] + 0.5 * forces[i][1] * timestep
            velocities[i][2] = half_velocities[i][2] + 0.5 * forces[i][2] * timestep

        if s % 10 == 0:
            potential_energy = calculate_potential_energy(positions, box_length, cut)
            temperature = calculate_temperature(velocities)
            kinetic_energy = 3 * n_particles * temperature
            total_energy = potential_energy + kinetic_energy
            total_energies.append(total_energy)
            print(f"Step: {s} \tTemperature: {temperature}\t Energy:{total_energy}")
            write_xyz(positions, f"trajectory_{n_particles}_{target_density}.xyz", "a")

    plt.plot(total_energies)
    plt.savefig("total_energies_.png")
