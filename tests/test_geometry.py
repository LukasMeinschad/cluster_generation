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



