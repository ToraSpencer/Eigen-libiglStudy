
// ������ת��Ϣ�õ���ת��������1����������ת����������ת�Ƕȣ�������ת���� 
template <typename ScalarVO, typename DerivedVA>
void getRotationMat(Eigen::Matrix<ScalarVO, 3, 3>& rotation, \
	const Eigen::PlainObjectBase<DerivedVA>& axisArrow, const double theta)
{
	assert((3 == axisArrow.rows() * axisArrow.cols()) && "assert!!! the input axisArrow should be in 3D space.");
	Eigen::Matrix<ScalarVO, 3, 1> axis = axisArrow.transpose().normalized().array().cast<ScalarVO>();
	rotation = Eigen::AngleAxis<ScalarVO>(theta, axis).toRotationMatrix();
}


// ������ת��Ϣ�õ���ת��������2�����õ���originArrow��ת��targetArrow����ת���� 
template <typename ScalarVO, typename DerivedV1, typename DerivedV2>
void getRotationMat(Eigen::Matrix<ScalarVO, 3, 3>& rotation, \
	const Eigen::PlainObjectBase<DerivedV1>& originArrow, \
	const Eigen::PlainObjectBase<DerivedV2>& targetArrow)
{
	assert((3 == originArrow.rows() * originArrow.cols()) && "assert!!! the input originArrow should be in 3D space.");
	assert((3 == targetArrow.rows() * targetArrow.cols()) && "assert!!! the input targetArrow should be in 3D space.");
	using RowVector3O = Eigen::Matrix<ScalarVO, 1, 3>;

	const double eps = 1e-5;
	rotation = Eigen::Matrix<ScalarVO, 3, 3>::Identity();
	if (0 == originArrow.norm() || 0 == targetArrow.norm())
		return;

	// 1. ��ԭ������Ŀ�귽��ƽ�У�
	RowVector3O originArrowCast = originArrow.array().cast<ScalarVO>();
	RowVector3O targetArrowCast = targetArrow.array().cast<ScalarVO>();
	RowVector3O axisArrowCast = originArrowCast.cross(targetArrowCast);			// ��ת�᣻
	if (std::abs(axisArrowCast.norm()) < eps)
	{
		// 1.1 ��ԭ������Ŀ�귽����ͬ��
		if (originArrowCast.dot(targetArrowCast) > 0) 
			return;
		else         // 1.2 ��ԭ������Ŀ�귽���෴��
		{
			// ����originArrowΪ�����ƽ���ڵ�һ��������
			axisArrowCast = RowVector3O{1, 1, 1};
			ScalarVO x0 = originArrowCast(0);
			ScalarVO y0 = originArrowCast(1);
			ScalarVO z0 = originArrowCast(2); 

			// ������o������Ϊx���������趨a������������������Ϊ1, x�������������o.dot(a)���x�������ٹ�һ���õ���ת������a
			if (std::abs(x0) >= std::abs(y0) && std::abs(x0) >= std::abs(z0))
				axisArrowCast(0) = -(y0 + z0)/x0;
			else if (std::abs(y0) >= std::abs(x0) && std::abs(y0) >= std::abs(z0))
				axisArrowCast(1) = -(x0 + z0) / y0;
			else
				axisArrowCast(2) = -(x0 + y0) / z0;
			axisArrowCast.normalize();
			ScalarVO x1 = axisArrowCast(0);
			ScalarVO y1 = axisArrowCast(1);
			ScalarVO z1 = axisArrowCast(2);

			rotation << 2*x1*x1 - 1, 2 * x1 * y1, 2 * x1*z1, \
						2*x1*y1, 2*y1*y1 -1, 2*y1*z1, \
						2*x1*z1, 2*y1*z1, 2*z1*z1 - 1;

			return;
		}
	}
		 

	axisArrowCast.normalize();
	ScalarVO x0 = axisArrowCast(0);
	ScalarVO y0 = axisArrowCast(1);
	ScalarVO z0 = axisArrowCast(2);
	ScalarVO cosTheta = originArrowCast.dot(targetArrowCast) / (originArrowCast.norm() * targetArrowCast.norm());
	ScalarVO sinTheta = std::sqrt(1 - cosTheta * cosTheta);

	// �ȼ���Eigen::AngleAxis<T>(theta, axis).toRotationMatrix()��������������������תtheta�Ƕȣ�
	rotation << cosTheta + (1 - cosTheta) * x0 * x0, (1 - cosTheta)* x0* y0 - sinTheta * z0, (1 - cosTheta)* x0* z0 + sinTheta * y0, \
		(1 - cosTheta)* y0* x0 + sinTheta * z0, cosTheta + (1 - cosTheta) * y0 * y0, (1 - cosTheta)* y0* z0 - sinTheta * x0, \
		(1 - cosTheta)* z0* x0 - sinTheta * y0, (1 - cosTheta)* z0* y0 + sinTheta * x0, cosTheta + (1 - cosTheta) * z0 * z0;
}