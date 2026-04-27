import numpy as np
import pytest

from modules.geometry import ReferenceFrame, Quaternion, GeometryOps

def test_reference_frame_construction():
    # Try constructing the dataclas
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([0.0, 1.0, 0.0])
    p3 = np.array([0.0, 0.0, 1.0])
    origin = np.array([0.0, 0.0, 0.0])
    frame = ReferenceFrame(axes=np.array([p1, p2, p3]), origin=origin)

    # Get the properties
    assert frame.x_axis is not None
    assert frame.y_axis is not None
    assert frame.z_axis is not None
    assert frame.origin is not None


@pytest.fixture
def initialize_quaternion_1():
    return Quaternion.from_axis_angle(np.array([0.0, 0.0, 1.0]), np.pi / 2)

@pytest.fixture
def initialize_quaternion_2():
    return Quaternion.from_axis_angle(np.array([0.0, 1.0, 0.0]), np.pi / 4)

def test_random_so3():
    for _ in range(10):
        R = Quaternion.random_so3()
        assert Quaternion.is_valid_rotation_matrix(R)

def test_from_two_vectors():
    """
    Test the function from_two_vectors that given two vectors,
    constructs a rotation quaternion to align the first vector to the second.
    """
    v1 = np.array([1.0, 0.0, 0.0])
    v2 = np.array([0.0, 1.0, 0.0])
    q = Quaternion.from_two_vectors(v1, v2)
    R = q.to_rotation_matrix()
    v1_rotated = R @ v1
    assert np.allclose(v1_rotated, v2)

def test_normalize(initialize_quaternion_1):
    """  
    Test the normalization of a quaternion. The norm should be 1 after normalization.
    """
    q = initialize_quaternion_1
    q.normalize()
    norm = np.linalg.norm(q.w + 1j*q.x + 1j*q.y + 1j*q.z)
    assert np.isclose(norm, 1.0)

def test_conjugate(initialize_quaternion_1):
    """  
    Tests the conjugate of a quaternion. The conjugate should have the same scalar part and negated vector part.
    """
    q = initialize_quaternion_1
    q_conj = q.conjugate()
    assert np.isclose(q.w, q_conj.w)  # Scalar part should be the same
    # Check x ,y, z (i,j,k) components are negated
    assert np.isclose(q.x, -q_conj.x)
    assert np.isclose(q.y, -q_conj.y)
    assert np.isclose(q.z, -q_conj.z)

def test_multiply(initialize_quaternion_1, initialize_quaternion_2):
    """  
    Tests the multiplication of two quaternions. The resulting quaternion should represent the combined rotation.
    """
    q1 = initialize_quaternion_1
    q2 = initialize_quaternion_2
    q_product = q1.multiply(q2)
    R1 = q1.to_rotation_matrix()
    R2 = q2.to_rotation_matrix()
    R_product = q_product.to_rotation_matrix()
    assert np.allclose(R_product, R1 @ R2)

def test_from_rotation_matrix(initialize_quaternion_1):
    """  
    Tests the conversion from a rotation matrix to a quaternion. The resulting quaternion should represent the same rotation as the original.
    """
    q = initialize_quaternion_1 
    R = q.to_rotation_matrix()
    q_from_R = Quaternion.from_rotation_matrix(R)
    R_from_q = q_from_R.to_rotation_matrix()
    assert np.allclose(R, R_from_q)

def test_rotate_point(initialize_quaternion_1):
    """  
    Tests the rotation of a point using a quaternion.
    For this the point is converted into a pure quaternion, then qpq^-1 is computed 
    """
    q = initialize_quaternion_1
    point = np.array([1.0, 0.0, 0.0])
    rotated_point = q.rotate_point(point)
    expected_rotated_point = np.array([0.0, 1.0, 0.0])  # 90 degree rotation around z-axis
    assert np.allclose(rotated_point, expected_rotated_point)

def test_rotate_points(initialize_quaternion_1):
    """  
    Tests the rotation of multiple points using a quaternion. The resulting points should be rotated correctly.
    """
    q = initialize_quaternion_1
    points = np.array([[1.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0],
                       [0.0, 0.0, 1.0]])
    rotated_points = q.rotate_points(points)
    expected_rotated_points = np.array([[0.0, 1.0, 0.0],
                                        [-1.0, 0.0, 0.0],
                                        [0.0, 0.0, 1.0]])  # 90 degree rotation around z-axis
    assert np.allclose(rotated_points, expected_rotated_points)

def test_centroid_compute():
    """  
    Tests the computation of a centroid using GeometryOps
    """
    points = np.array([[1.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0],
                       [0.0, 0.0, 1.0]])
    centroid = GeometryOps.centroid(points)
    expected_centroid = np.array([1/3, 1/3, 1/3])
    assert np.allclose(centroid, expected_centroid)

def test_center_of_mass_compute():
    """  
    Tests the computation of a center of mass using GeometryOps
    """
    points = np.array([[1.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0],
                       [0.0, 0.0, 1.0]])
    masses = np.array([1.0, 2.0, 3.0])
    com = GeometryOps.center_of_mass(points, masses)
    expected_com = np.array([1/6, 1/3, 1/2])
    assert np.allclose(com, expected_com)

def test_inertia_tensor_compute():
    """  
    Tests the computation of the moment of inertia tensor using GeometryOp
    """
    points = np.array([[1.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0],
                       [0.0, 0.0, 1.0]])
    masses = np.array([1.0, 1.0, 1.0])
    inertia_tensor = GeometryOps.inertia_tensor(points, masses)
    # The inertia tensor would be a matrix with 2 on the diagonal
    expected_inertia_tensor = np.array([[2.0, 0.0, 0.0],
                                        [0.0, 2.0, 0.0],
                                        [0.0, 0.0, 2.0]])
    assert np.allclose(inertia_tensor, expected_inertia_tensor)

def test_kabsch_rotation():
    """ 
    Tests the Kabsch algorithm for finding the optimal rotation between point P,Q
    """ 
    P = np.array([[1.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0],
                  [0.0, 0.0, 1.0]])
    Q = np.array([[0.0, 1.0, 0.0],
                  [-1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0]])  # P rotated by 90 degrees around z-axis
    R = GeometryOps.kabsch_rotation(P, Q)
    # Kabsch returns P @ R = Q 
    P_rotated = P @ R
    assert np.allclose(P_rotated, Q)

def test_find_optimal_correspondence():
    """  
    Tests the optimal correspondence function between two sets of coordinates P and Q using the
    Hungarian algorithm to solve the linear assignment problem
    """
    P = np.array([[1.0,2.0,3.0], [4.0,5.0,6.0], [7.0,8.0,9.0]])
    Q = P.copy()
    result = GeometryOps.find_optimal_correspondence(P, Q)
    # Since P and Q are the same, the optimal correspondence should be the identity mapping
    assert np.allclose(result, P)

def test_pairwise_distance_matrix():
    """  
    Test pairwise distance function 
    """
    # Computes the pairwise distance matrix for a set of points
    points = np.array([[0.0, 0.0, 0.0],
                       [1.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0]])
    distance_matrix = GeometryOps.pairwise_distance_matrix(points)
    expected_distance_matrix = np.array([[0.0, 1.0, 1.0],
                                         [1.0, 0.0, np.sqrt(2)],
                                         [1.0, np.sqrt(2), 0.0]])
    assert np.allclose(distance_matrix, expected_distance_matrix)

def test_normalize_vector():
    """  
    Test the normalization function of GeometryOps
    """
    vector = np.array([1.0, 1.0, 1.0])
    normalized_vector = GeometryOps.normalize_vector(vector)
    expected_normalized_vector = np.array([1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)])
    assert np.allclose(normalized_vector, expected_normalized_vector)

def test_rotation_matrix_rodrigues():
    """   
    Tests the construction of a rotation matrix using Rodrigues formula
    """
    axis = np.array([0.0, 0.0, 1.0])
    angle = np.pi / 2
    R = GeometryOps.rotation_matrix_rodrigues(axis, angle)
    expected_R = np.array([[0.0, -1.0, 0.0],
                           [1.0, 0.0, 0.0],
                           [0.0, 0.0, 1.0]])  # 90 degree rotation around z-axis
    assert np.allclose(R, expected_R)

def test_compute_optimal_correspondence_rmsd():
    """  
    Tests the computation of the optimal correspondence RMSD between two sets of coordinates P and Q
    """
    P = np.array([[1.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0],
                  [0.0, 0.0, 1.0]])
    Q = np.array([[0.0, 1.0, 0.0],
                  [-1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0]])  # P rotated by 90 degrees around z-axis
    rmsd = GeometryOps.compute_optimal_correspondence_rmsd(P, Q)
    assert np.isclose(rmsd, 0.0)

def test_align_vectors():
    """  
    Tests the alignment function using the Rodrigues formula for the rotation
    """
    v = np.array([1.0, 0.0, 0.0])
    target = np.array([0.0, 1.0, 0.0])
    R = GeometryOps.align_vectors(v, target)
    v_aligned = R @ v
    assert np.allclose(v_aligned, target)

def test_distance():
    """  
    Tests the distance function between two points
    """
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([0.0, 1.0, 0.0])
    distance = GeometryOps.distance(p1, p2)
    expected_distance = np.sqrt(2)
    assert np.isclose(distance, expected_distance)




