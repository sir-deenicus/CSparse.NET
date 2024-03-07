/*
  
         
*/
namespace CSparse.Single
{
    using System;

    /// <summary>
    /// Vector helper methods.
    /// </summary>
    
    public static class Vector
    { 
        /// <summary>
        /// Copy one vector to another.
        /// </summary>
        /// <param name="src">The source array.</param>
        /// <param name="dst">The destination array.</param>
        public static void Copy(float[] src, float[] dst)
        {
            Buffer.BlockCopy(src, 0, dst, 0, src.Length * Constants.SizeOfSingle);
        }

        /// <summary>
        /// Copy one vector to another.
        /// </summary>
        /// <param name="n">Number of values to copy.</param>
        /// <param name="src">The source array.</param>
        /// <param name="dst">The destination array.</param>
        public static void Copy(int n, float[] src, float[] dst)
        {
            Buffer.BlockCopy(src, 0, dst, 0, n * Constants.SizeOfSingle);
        }

        /// <summary>
        /// Create a new vector.
        /// </summary>
        public static float[] Create(int length, float value)
        {
            float[] result = new float[length];

            for (int i = 0; i < length; i++)
            {
                result[i] = value;
            }

            return result;
        }

        /// <summary>
        /// Clone the given vector.
        /// </summary>
        public static float[] Clone(float[] src)
        {
            float[] result = new float[src.Length];

            Buffer.BlockCopy(src, 0, result, 0, src.Length * Constants.SizeOfSingle);

            return result;
        }

        /// <summary>
        /// Set vector values to zero.
        /// </summary>
        public static void Clear(float[] x)
        {
            Array.Clear(x, 0, x.Length);
        }

        /// <summary>
        /// Computes the dot product of two vectors.
        /// </summary>
        [Obsolete("Use DotProduct(n, x, y).")]
        public static float DotProduct(float[] x, float[] y)
        {
            return DotProduct(x.Length, x, y);
        }

        /// <summary>
        /// Computes the dot product of two vectors.
        /// </summary>
        public static float DotProduct(int n, float[] x, float[] y)
        {
            float result = 0.0f;

            for (int i = 0; i < n; i++)
            {
                result += x[i] * y[i];
            }

            return result;
        }

        /// <summary>
        /// Computes the pointwise product of two vectors.
        /// </summary>
        [Obsolete("Use PointwiseMultiply(n, x, y, target).")]
        public static void PointwiseMultiply(float[] x, float[] y, float[] target)
        {
            PointwiseMultiply(x.Length, x, y, target);
        }

        /// <summary>
        /// Computes the pointwise product of two vectors.
        /// </summary>
        public static void PointwiseMultiply(int n, float[] x, float[] y, float[] target)
        {
            for (int i = 0; i < n; i++)
            {
                target[i] = x[i] * y[i];
            }
        }

        /// <summary>
        /// Computes the norm of a vector, sqrt( x' * x ).
        /// </summary>
        [Obsolete("Use Norm(n, x).")]
        public static float Norm(float[] x)
        {
            return Norm(x.Length, x);
        }

        /// <summary>
        /// Computes the norm of a vector, sqrt( x' * x ).
        /// </summary>
        public static float Norm(int n, float[] x)
        {
            float result = 0.0f;

            for (int i = 0; i < n; ++i)
            {
                result += x[i] * x[i];
            }

            return (float)Math.Sqrt(result);
        }

        /// <summary>
        /// Computes the norm of a vector avoiding overflow, sqrt( x' * x ).
        /// </summary>
        [Obsolete("Use NormRobust(n, x).")]
        public static float NormRobust(float[] x)
        {
            return NormRobust(x.Length, x);
        }

        /// <summary>
        /// Computes the norm of a vector avoiding overflow, sqrt( x' * x ).
        /// </summary>
        public static float NormRobust(int n, float[] x)
        {
            float scale = 0.0f, ssq = 1.0f;

            for (int i = 0; i < n; ++i)
            {
                if (x[i] != 0.0f)
                {
                    float absxi = Math.Abs(x[i]);
                    if (scale < absxi)
                    {
                        ssq = 1.0f + ssq * (scale / absxi) * (scale / absxi);
                        scale = absxi;
                    }
                    else
                    {
                        ssq += (absxi / scale) * (absxi / scale);
                    }
                }
            }

            return scale * (float)Math.Sqrt(ssq);
        }

        /// <summary>
        /// Scales a vector by a given factor, x = a * x.
        /// </summary>
        public static void Scale(float a, float[] x)
        {
            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                x[i] *= a;
            }
        }

        /// <summary>
        /// Scales a vector by a given factor, target = a * x.
        /// </summary>
        public static void Scale(int n, float a, float[] x, float[] target)
        {
            for (int i = 0; i < n; i++)
            {
                target[i] = a * x[i];
            }
        }

        /// <summary>
        /// Add a scaled vector to another vector, y = a * x + y.
        /// </summary>
        public static void Axpy(float a, float[] x, float[] y)
        {
            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                y[i] += a * x[i];
            }
        }

        /// <summary>
        /// Add two vectors, z = a * x + b * y.
        /// </summary>
        public static void Add(float a, float[] x, float b, float[] y, float[] target)
        {
            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                target[i] = a * x[i] + b * y[i];
            }
        }

        /// <summary>
        /// Add two vectors, target = a * x + y.
        /// </summary>
        public static void Add(int n, float a, float[] x, float[] y, float[] target)
        {
            for (int i = 0; i < n; i++)
            {
                target[i] = a * x[i] + y[i];
            }
        }

        /// <summary>
        /// Add two vectors, target = a * x + b * y.
        /// </summary>
        public static void Add(int n, float a, float[] x, float b, float[] y, float[] target)
        {
            for (int i = 0; i < n; i++)
            {
                target[i] = a * x[i] + b * y[i];
            }
        }

        /// <summary>
        /// Add three vectors, target = a * x + b * y + z.
        /// </summary>
        public static void Add(int n, float a, float[] x, float b, float[] y, float[] z, float[] target)
        {
            for (int i = 0; i < n; i++)
            {
                target[i] = a * x[i] + b * y[i] + z[i];
            }
        }

        /// <summary>
        /// Add three vectors, target = a * x + b * y + c * z.
        /// </summary>
        public static void Add(int n, float a, float[] x, float b, float[] y, float c, float[] z, float[] target)
        {
            for (int i = 0; i < n; i++)
            {
                target[i] = a * x[i] + b * y[i] + c * z[i];
            }
        }

    }

}