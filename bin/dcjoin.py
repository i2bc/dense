import collections
import itertools


class AbstractMarker:
    """
    General Representation for a marker in DCJ model.
    """
    # TODO Should look into abstract/metaclasses if appropriate
    def __init__(self):
        pass


class Telomere(AbstractMarker):
    """
    Representation for a telomere in DCJ model.
    """
    def __init__(self):
        self.reverse = False
        pass

    def __repr__(self):
        return "."

    def __eq__(self, other):
        return type(self) == type(other)

    def __hash__(self):
        return hash(".")


class Marker(AbstractMarker):
    """
    Representation for DNA fragment in DCJ model.

    Parameters
    ----------
    uid: int (required)
        unique id

    dna: str
        What dna fragment does this marker represent?

    reverse: boolean
        Is the DNA fragment reversed?
    """

    def __init__(self, uid, dna="", reverse=False):
        self.reverse = reverse
        self.dna = dna
        self.uid = uid

    def __eq__(self, other):
        if type(other) is type(self):
            return self.uid == other.uid
        return False

    def __repr__(self):
        return "-" + str(self.uid) if self.reverse else str(self.uid)

    def __hash__(self):
        return hash(self.uid)


class MarkerEnd:
    """Markers with designated head or tail end. Used in adjacencies"""
    def __init__(self, marker, head=True):
        self.marker = marker
        self.head = head

    def __repr__(self):
        return str(self.marker) if self.head else str(self.marker) + "*"

    def __eq__(self, other):
        return self.marker == other.marker and self.head == other.head

    def other_adjacency_end(self, adj):
        other_end_marker = adj.right_end_marker if adj.left_end_marker == self else adj.left_end_marker
        return other_end_marker

    def __hash__(self):
        return hash(self.marker) * hash(self.head)


class Chromosome(collections.UserList):
    """The abstraction of a chromosome in DCJ model, a list of telomeres and markers

    Parameters
    ----------
    markers: list or string (req)
        The object to be converted into a chromosome.

    """
    # Just a list of markers

    string_d = {}
    uid_counter = 0

    def __init__(self, markers):
        self.marker_set = set()
        # If it's already a list of markers, skip all this. Otherwise.
        if isinstance(markers, list):
            for marker in markers:
                if not isinstance(marker, AbstractMarker):
                    raise TypeError("Chromosome must be a list of Marker instances or a string of unique markers.")
            content = markers
        else:
        # Otherwise loop through like it is a string
            content = []
            # If first is telomere, last must be telomere or vice-versa
            if not (markers[0] == ".") == (markers[-1] == "."):
                raise ValueError("Linear Chromosome must start and end with telomeres.")
            for i, s in enumerate(markers):
                dna = s
                rev = s.isupper()
                telomere = s == "."
                if s.lower() in self.marker_set and not telomere:
                    raise ValueError("Duplicated markers are not allowed. ({0})".format(s.lower()))
                elif s.lower() in Chromosome.string_d:
                    uid = Chromosome.string_d[s.lower()]
                    self.marker_set.add(s.lower())
                else:
                    Chromosome.string_d[s.lower()] = Chromosome.uid_counter
                    self.marker_set.add(s.lower())
                    uid = Chromosome.uid_counter
                    Chromosome.uid_counter += 1
                if telomere and i not in [0, len(markers) - 1]:
                    raise ValueError("Telomere cannot appear in middle of chromosome.")
                elif telomere:
                    content.append(Telomere())
                else:
                    content.append(Marker(uid, dna, rev))
        super().__init__(content)

    @classmethod
    def clear(cls):
        """
        Resets the string_d and uid_counter.
        :return:
        """

        cls.string_d = {}
        cls.uid_counter = 0


class Genome(collections.UserList):
    """The abstraction of a genome in DCJ model, a list of chromosomes"""
    def __init__(self, *args):
        for c in args:
            if not isinstance(c, Chromosome):
                raise TypeError("Argument list passed to Genome constructor must be Chromosome class instance")

        if len(args) > 1:
            duplicated_markers = set.intersection(*[c.marker_set for c in args])
        else:
            duplicated_markers = None

        if duplicated_markers:
            raise ValueError("Duplicated markers found across chromosomes. ({0})".format(duplicated_markers))
        super().__init__(args)


class Adjacency(collections.UserList):
    """Adjacency Data Structure"""
    def __init__(self, left_end_marker=None, right_end_marker=None, label=None):
        # List of markers, chromosome.
        self.label = label
        self.left_end_marker = left_end_marker
        self.right_end_marker = right_end_marker
        super().__init__([left_end_marker, right_end_marker])

    def __repr__(self):
        return "[" + str(self.left_end_marker) + ", " + str(self.label) + ", " + str(self.right_end_marker) + "]"


class AdjacencyGraph:
    """Adjacency information for a genome"""
    # General Idea
    # Find common elements
    # Make Adjacency lists
    # Make two tables for fast indexing
    # Create Adjacency Graph
    # Traverse Graph counting cycles, paths and runs

    def __init__(self, A, B):
        # PseudoCode
        # Find A, B, G marker set
        # Create Adjacency lists

        # Finding Common Elements
        marker_set_a = set(itertools.chain.from_iterable(A))
        marker_set_b = set(itertools.chain.from_iterable(B))
        self.commonMarkers = marker_set_a.intersection(marker_set_b)
        if Telomere() in self.commonMarkers:
            self.commonMarkers.remove(Telomere())
        marker_set_a = marker_set_a.difference(self.commonMarkers)
        marker_set_b = marker_set_b.difference(self.commonMarkers)

        # Make Adjacency lists
        self.adjA = []
        self.adjB = []
        adj = [None, None]
        # Label is a list of markers NOT endmarkers

        # Created Adjacencies for Genomes A and B
        adjacencies = [self.adjA, self.adjB]
        reference_A = {MarkerEnd(Telomere(), head=True): None}
        reference_B = {MarkerEnd(Telomere(), head=True): None}
        references = [reference_A, reference_B]
        for i, genome in enumerate([A, B]):
            adjacency_length = 0
            for chromosome in genome:
                index = 0
                adj_index = 0
                current_marker = chromosome[index]
                chromosome_markers = set([x for x in chromosome])
                adjacency_length += len(self.commonMarkers.intersection(chromosome_markers))
                if chromosome[0] == Telomere():
                    adjacency_length += 1
                while len(adjacencies[i]) < adjacency_length:
                    label = []

                    # Fill first adjacency position
                    if current_marker == Telomere():
                        adj[0] = MarkerEnd(Telomere(), head=True)
                    else:
                        adj[0] = MarkerEnd(current_marker, head=not current_marker.reverse)

                    # While loop to add to label
                    next_marker = chromosome[0] if index >= len(chromosome) - 1 else chromosome[index+1]
                    while next_marker not in self.commonMarkers:
                        # Fixes linear chromosomes
                        if next_marker == Telomere():
                            break
                        label.append(next_marker)
                        index += 1
                        # Circular genome, must wrap around
                        next_marker = chromosome[0] if index >= len(chromosome) - 1 else chromosome[index + 1]

                    # fill second adjacency position
                    if next_marker == Telomere():
                        adj[1] = MarkerEnd(Telomere(), head=True)
                    else:
                        adj[1] = MarkerEnd(next_marker, head=next_marker.reverse)
                    adjacencies[i].append(Adjacency(left_end_marker=adj[0], right_end_marker=adj[1], label=label))

                    # Create reference matrix for A and B
                    if adj[0].marker != Telomere():
                        references[i][adj[0]] = len(adjacencies[i]) - 1
                    if adj[1].marker != Telomere():
                        references[i][adj[1]] = len(adjacencies[i]) - 1

                    current_marker = next_marker
                    index += 1
                    adj_index += 1

        # Traverse the Adjacency graph
        to_visit_a_index = set(range(len(self.adjA)))
        to_visit_b_index = set(range(len(self.adjB)))
        visited_a_index = set()
        visited_b_index = set()

        cycles = 0
        ab_paths = 0
        a_runs = 0
        on_a_run = False
        b_runs = 0
        on_b_run = False

        while to_visit_a_index:
            current_adj_index = to_visit_a_index.pop()
            visited_a_index.add(current_adj_index)

            left_marker = self.adjA[current_adj_index].left_end_marker
            right_marker = self.adjA[current_adj_index].right_end_marker

            # Check label for A_run
            if self.adjA[current_adj_index].label:
                a_runs += 1
                on_a_run = True

            # Checks for AB paths
            paths_end_on_a = [True, True]
            for i, current_marker in enumerate([left_marker, right_marker]):
                a_side = True
                next_adj_index = reference_B[current_marker]
                if next_adj_index in visited_b_index:
                    # If cycle detected on left_marker, no need to check right marker
                    continue
                while next_adj_index is not None:
                    current_adj_index = next_adj_index
                    a_side = not a_side
                    adj_side = self.adjA if a_side else self.adjB
                    current_adj = adj_side[current_adj_index]
                    reference_side = reference_B if a_side else reference_A

                    visited_a_index.add(current_adj_index) if a_side else visited_b_index.add(current_adj_index)
                    to_visit_a_index.remove(current_adj_index) if a_side else to_visit_b_index.remove(current_adj_index)

                    # Check label, and update run logs
                    if current_adj.label:
                        if not a_side and on_a_run:
                            on_b_run = True
                            on_a_run = False
                            b_runs += 1
                        elif not a_side and not on_b_run:
                            on_b_run = True
                            b_runs += 1

                        elif a_side and on_b_run:
                            on_a_run = True
                            on_b_run = False
                            a_runs += 1
                        elif a_side and not on_a_run:
                            on_a_run = True
                            a_runs += 1


                    # Switch markers on adjacency
                    if current_marker != current_adj.left_end_marker:
                        current_marker = current_adj.left_end_marker
                    else:
                        current_marker = current_adj.right_end_marker

                    # Find next pointer location
                    next_adj_index = reference_side[current_marker]

                    if next_adj_index is None:
                        # Log which genome path ends on
                        paths_end_on_a[i] = a_side
                        break
                    elif (a_side and next_adj_index in visited_b_index)\
                            or (not a_side and next_adj_index in visited_a_index):
                        cycles += 1
                        break
                    else:
                        # The path continues
                        pass

            # AB_path is detected
            if paths_end_on_a[0] != paths_end_on_a[1]:
                ab_paths += 1

        self.cycles = cycles
        self.ab_paths = ab_paths
        self.a_runs = a_runs
        self.b_runs = b_runs
        self.run_potential = a_runs + b_runs
        self.indel_potential = (self.run_potential + 1) // 2 + ((self.run_potential // 2) % 2) if self.run_potential > 0 else 0


def calculate_distance(A, B, method="dcj"):
    """Calculate the DCJ distance between genome A and B"""
    if method not in ["dcj", "indel"]:
        raise ValueError('method argument must be "dcj" (default) or "indel"')

    # DCJ distance
    # DCJ(A, B) = |G| - c - b/2
    # clear dictionary for chromosomes
    Chromosome.clear()
    genome_a = transform_genome(A)
    genome_b = transform_genome(B)

    ag = AdjacencyGraph(genome_a, genome_b)
    dcj_distance = len(ag.commonMarkers) - ag.cycles - ag.ab_paths / 2
    dcj_indel_distance = dcj_distance + ag.indel_potential

    if method == "indel":
        return dcj_indel_distance
    return dcj_distance


def transform_genome(genome):
    """
    Transforms the genome into the underlying class structure.

    :param genome:
    :return:
    """
    return genome if isinstance(genome, Genome) else Genome(*[Chromosome(c) for c in genome])

if __name__ == '__main__':
    print(calculate_distance([".abc."], [".aCBd."], method="dcj"))
    pass

