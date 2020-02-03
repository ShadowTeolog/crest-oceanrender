using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Numerics;
using WaveShapeCalculator.Random;

namespace CrestServer
{
    internal static class Mathf
    {
        public const float Deg2Rad=(float) (Math.PI/180.0);
        public const float PI=(float) Math.PI;
        public static float Pow(float f, float f1) => (float) Math.Pow(f, f1);
        public static float Log(float wavelength) => (float) Math.Log(wavelength);
        public static float Clamp(float value, float min, float max) => Math.Min(max, Math.Max(min, value));
        public static float Sqrt(float a2) => (float) Math.Sqrt(a2);
        public static float Floor(float wlPow2) => (float) Math.Floor(wlPow2);
        public static float Log10(float pow)=>(float) Math.Log10(pow);
        public static float Exp(float f) => (float)Math.Exp(f);

        public static float Lerp(float start, float end, float ratio)
        {
            ratio = Clamp(ratio, 0, 1);
            return start * (1 - ratio) + end * ratio;
        }

        public static float Cos(float angleRadians) => (float) Math.Cos(angleRadians);
        public static float Sin(float angleRadians) => (float) Math.Sin(angleRadians);

        public static float Repeat(float f, float f1) => f % f1;
    }


    /// <summary>
    /// Support script for Gerstner wave ocean shapes.
    /// Generates a number of batches of Gerstner waves.
    /// </summary>
    public class ShapeGerstnerBatched :  ICollProvider
    {
        
        //[Tooltip("The spectrum that defines the ocean surface shape. Create asset of type Crest/Ocean Waves Spectrum.")]
        public OceanWaveSpectrum _spectrum;

        //[Delayed, Tooltip("How many wave components to generate in each octave.")]
        public int _componentsPerOctave = 5;

        //[Range(0f, 1f)]
        public float _weight = 1f;
        public uint _randomSeed = 0;
        public float windAngle;

        // data for all components
        float[] _wavelengths;
        float[] _amplitudes;
        float[] _angleDegs;
        float[] _phases;

        private MersenneTwister randomGenerator;
        // useful references
        //Material[] _materials;
        bool[] _drawLOD;
        //Material _materialBigWaveTransition;
        bool _drawLODTransitionWaves;

        // Shader to be used to render evaluate Gerstner waves for each LOD
        //Shader _waveShader;

        // IMPORTANT - this mirrors the constant with the same name in ShapeGerstnerBatch.shader, both must be updated together!
        const int BATCH_SIZE = 32;

        enum CmdBufStatus
        {
            NoStatus,
            NotAttached,
            Attached
        }

        // scratch data used by batching code
        struct UpdateBatchScratchData
        {
            public static Vector4[] _twoPiOverWavelengthsBatch = new Vector4[BATCH_SIZE / 4];
            public static Vector4[] _ampsBatch = new Vector4[BATCH_SIZE / 4];
            public static Vector4[] _waveDirXBatch = new Vector4[BATCH_SIZE / 4];
            public static Vector4[] _waveDirZBatch = new Vector4[BATCH_SIZE / 4];
            public static Vector4[] _phasesBatch = new Vector4[BATCH_SIZE / 4];
            public static Vector4[] _chopAmpsBatch = new Vector4[BATCH_SIZE / 4];
        }

        public ShapeGerstnerBatched()
        {
            _spectrum = new OceanWaveSpectrum();
        }

        
        void InitPhases(System.Random random)
        {
            // Set random seed to get repeatable results
            
            var totalComps = _componentsPerOctave * OceanWaveSpectrum.NUM_OCTAVES;
            _phases = new float[totalComps];
            for (var octave = 0; octave < OceanWaveSpectrum.NUM_OCTAVES; octave++)
            {
                for (var i = 0; i < _componentsPerOctave; i++)
                {
                    var index = octave * _componentsPerOctave + i;
                    var rnd = (i + (float)random.NextDouble()) / _componentsPerOctave;
                    _phases[index] = 2f * Mathf.PI * rnd;
                }
            }
        }

        public void SetOrigin(Vector3 newOrigin)
        {
            if (_phases == null) return;
            for (int i = 0; i < _phases.Length; i++)
            {
                var direction = new Vector3(Mathf.Cos((windAngle + _angleDegs[i]) * Mathf.Deg2Rad), 0f, Mathf.Sin((windAngle + _angleDegs[i]) * Mathf.Deg2Rad));
                var phaseOffsetMeters = Vector3.Dot(newOrigin, direction);
                // wave number
                var k = 2f * Mathf.PI / _wavelengths[i];

                _phases[i] = Mathf.Repeat(_phases[i] + phaseOffsetMeters * k, Mathf.PI * 2f);
            }
        }

        public void Update()
        {
            //create or reseed
            if (randomGenerator == null)
                randomGenerator = new MersenneTwister(_randomSeed);
            else
                randomGenerator.Seed = _randomSeed;
            if (_phases == null || _phases.Length != _componentsPerOctave * OceanWaveSpectrum.NUM_OCTAVES)
            {
                InitPhases(randomGenerator);
                randomGenerator.Seed = _randomSeed; //reseed again after phases init
            }
            
            
            _spectrum.GenerateWaveData(randomGenerator,_componentsPerOctave, ref _wavelengths, ref _angleDegs);
            UpdateAmplitudes();
        }

        void UpdateAmplitudes()
        {
            if (_amplitudes == null || _amplitudes.Length != _wavelengths.Length)
            {
                _amplitudes = new float[_wavelengths.Length];
            }

            for (int i = 0; i < _wavelengths.Length; i++)
            {
                _amplitudes[i] = _weight * _spectrum.GetAmplitude(_wavelengths[i], _componentsPerOctave);
            }
        }

        
        float ComputeWaveSpeed(float wavelength/*, float depth*/)
        {
            // wave speed of deep sea ocean waves: https://en.wikipedia.org/wiki/Wind_wave
            // https://en.wikipedia.org/wiki/Dispersion_(water_waves)#Wave_propagation_and_dispersion
            float g = 9.81f;
            float k = 2f * Mathf.PI / wavelength;
            //float h = max(depth, 0.01);
            //float cp = sqrt(abs(tanh_clamped(h * k)) * g / k);
            float cp = Mathf.Sqrt(g / k);
            return cp;
        }

        public bool GetSurfaceVelocity(ref Vector3 i_worldPos, SamplingData i_samplingData,float time, out Vector3 o_surfaceVel)
        {
            o_surfaceVel = Vector3.Zero;

            if (_amplitudes == null) return false;

            Vector2 pos = new Vector2(i_worldPos.X, i_worldPos.Z);
            
            
            float minWaveLength = i_samplingData._minSpatialLength / 2f;

            for (int j = 0; j < _amplitudes.Length; j++)
            {
                if (_amplitudes[j] <= 0.001f) continue;
                if (_wavelengths[j] < minWaveLength) continue;

                float C = ComputeWaveSpeed(_wavelengths[j]);

                // direction
                Vector2 D = new Vector2(Mathf.Cos((windAngle + _angleDegs[j]) * Mathf.Deg2Rad), Mathf.Sin((windAngle + _angleDegs[j]) * Mathf.Deg2Rad));
                // wave number
                float k = 2f * Mathf.PI / _wavelengths[j];

                float x = Vector2.Dot(D, pos);
                float t = k * (x + C * time) + _phases[j];
                float disp = -_spectrum._chop * k * C * Mathf.Cos(t);
                o_surfaceVel += _amplitudes[j] * new Vector3(
                    D.X * disp,
                    -k * C * Mathf.Sin(t),
                    D.Y * disp
                    );
            }

            return true;
        }

        public bool SampleHeight(ref Vector3 i_worldPos, SamplingData i_samplingData,float time, out float o_height)
        {
            o_height = 0f;

            Vector3 posFlatland = i_worldPos;

            Vector3 undisplacedPos;
            if (!ComputeUndisplacedPosition(ref posFlatland, i_samplingData,time, out undisplacedPos))
                return false;

            Vector3 disp;
            if (!SampleDisplacement(ref undisplacedPos, i_samplingData,time, out disp))
                return false;

            o_height = posFlatland.Y + disp.Y;

            return true;
        }

        public void ReturnSamplingData(SamplingData i_data)
        {
            i_data._minSpatialLength = -1f;
        }

        public bool ComputeUndisplacedPosition(ref Vector3 i_worldPos, SamplingData i_samplingData,float time, out Vector3 o_undisplacedWorldPos)
        {
            // FPI - guess should converge to location that displaces to the target position
            Vector3 guess = i_worldPos;
            // 2 iterations was enough to get very close when chop = 1, added 2 more which should be
            // sufficient for most applications. for high chop values or really stormy conditions there may
            // be some error here. one could also terminate iteration based on the size of the error, this is
            // worth trying but is left as future work for now.
            Vector3 disp;
            for (int i = 0; i < 4 && SampleDisplacement(ref guess, i_samplingData,time, out disp); i++)
            {
                Vector3 error = guess + disp - i_worldPos;
                guess.X -= error.X;
                guess.Z -= error.Z;
            }

            o_undisplacedWorldPos = guess;
            return true;
        }

        public AvailabilityResult CheckAvailability(ref Vector3 i_worldPos, SamplingData i_samplingData)
        {
            return _amplitudes == null ? AvailabilityResult.NotInitialisedYet : AvailabilityResult.DataAvailable;
        }

        // Compute normal to a surface with a parameterization - equation 14 here: http://mathworld.wolfram.com/NormalVector.html
        public bool SampleNormal(ref Vector3 i_undisplacedWorldPos, SamplingData i_samplingData,float time, out Vector3 o_normal)
        {
            o_normal = Vector3.Zero;

            if (_amplitudes == null) return false;

            var pos = new Vector2(i_undisplacedWorldPos.X, i_undisplacedWorldPos.Z);
            float minWaveLength = i_samplingData._minSpatialLength / 2f;

            // base rate of change of our displacement function in x and z is unit
            var delfdelx = new Vector3(1, 0, 0);
            var delfdelz = new Vector3(0, 0, 1);

            for (int j = 0; j < _amplitudes.Length; j++)
            {
                if (_amplitudes[j] <= 0.001f) continue;
                if (_wavelengths[j] < minWaveLength) continue;

                float C = ComputeWaveSpeed(_wavelengths[j]);

                // direction
                var D = new Vector2(Mathf.Cos((windAngle + _angleDegs[j]) * Mathf.Deg2Rad), Mathf.Sin((windAngle + _angleDegs[j]) * Mathf.Deg2Rad));
                // wave number
                float k = 2f * Mathf.PI / _wavelengths[j];

                float x = Vector2.Dot(D, pos);
                float t = k * (x + C * time) + _phases[j];
                float disp = k * -_spectrum._chop * Mathf.Cos(t);
                float dispx = D.X * disp;
                float dispz = D.Y * disp;
                float dispy = -k * Mathf.Sin(t);

                delfdelx += _amplitudes[j] * new Vector3(D.X * dispx, D.X * dispy, D.Y * dispx);
                delfdelz += _amplitudes[j] * new Vector3(D.X * dispz, D.Y * dispy, D.Y * dispz);
            }

            o_normal = Vector3.Normalize(Vector3.Cross(delfdelz, delfdelx));

            return true;
        }

        public bool SampleDisplacement(ref Vector3 i_worldPos, SamplingData i_samplingData,float time, out Vector3 o_displacement)
        {
            o_displacement = Vector3.Zero;

            if (_amplitudes == null)
            {
                return false;
            }

            Vector2 pos = new Vector2(i_worldPos.X, i_worldPos.Z);
            float minWavelength = i_samplingData._minSpatialLength / 2f;

            for (int j = 0; j < _amplitudes.Length; j++)
            {
                if (_amplitudes[j] <= 0.001f) continue;
                if (_wavelengths[j] < minWavelength) continue;

                float C = ComputeWaveSpeed(_wavelengths[j]);

                // direction
                Vector2 D = new Vector2(Mathf.Cos((windAngle + _angleDegs[j]) * Mathf.Deg2Rad), Mathf.Sin((windAngle + _angleDegs[j]) * Mathf.Deg2Rad));
                // wave number
                float k = 2f * Mathf.PI / _wavelengths[j];

                float x = Vector2.Dot(D, pos);
                float t = k * (x + C * time) + _phases[j];
                float disp = -_spectrum._chop * Mathf.Sin(t);
                o_displacement += _amplitudes[j] * new Vector3(
                    D.X * disp,
                    Mathf.Cos(t),
                    D.Y * disp
                    );
            }

            return true;
        }

        public void SampleDisplacementVel(ref Vector3 i_worldPos, SamplingData i_samplingData,float time, out Vector3 o_displacement, out bool o_displacementValid, out Vector3 o_displacementVel, out bool o_velValid)
        {
            o_displacementValid = SampleDisplacement(ref i_worldPos, i_samplingData,time, out o_displacement);
            o_velValid = GetSurfaceVelocity(ref i_worldPos, i_samplingData,time, out o_displacementVel);
        }
    }
}

