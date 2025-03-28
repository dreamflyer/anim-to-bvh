(function webpackUniversalModuleDefinition(root, factory) {
	if(typeof exports === 'object' && typeof module === 'object')
		module.exports = factory();
	else if(typeof define === 'function' && define.amd)
		define([], factory);
	else if(typeof exports === 'object')
		exports["AnimToBvh"] = factory();
	else
		root["AnimToBvh"] = factory();
})(typeof self !== 'undefined' ? self : this, () => {
return /******/ (() => { // webpackBootstrap
/******/ 	"use strict";
/******/ 	var __webpack_modules__ = ({

/***/ "./node_modules/quaternion/dist/quaternion.js":
/*!****************************************************!*\
  !*** ./node_modules/quaternion/dist/quaternion.js ***!
  \****************************************************/
/***/ ((module) => {



/**
 * Creates a new Quaternion object
 * 
 * @param {number} w 
 * @param {number} x 
 * @param {number} y 
 * @param {number} z 
 * @returns 
 */
function newQuaternion(w, x, y, z) {
  const f = Object.create(Quaternion.prototype);

  f['w'] = w;
  f['x'] = x;
  f['y'] = y;
  f['z'] = z;

  return f;
}

/**
 * Creates a new normalized Quaternion object
 *
 * @param {number} w
 * @param {number} x
 * @param {number} y
 * @param {number} z
 * @returns
 */
function newNormalized(w, x, y, z) {
  const f = Object.create(Quaternion.prototype);

  // We assume |Q| > 0 for internal usage
  const il = 1 / Math.sqrt(w * w + x * x + y * y + z * z);

  f['w'] = w * il;
  f['x'] = x * il;
  f['y'] = y * il;
  f['z'] = z * il;

  return f;
}

/**
 * Calculates log(sqrt(a^2+b^2)) in a way to avoid overflows
 *
 * @param {number} a
 * @param {number} b
 * @returns {number}
 */
function logHypot(a, b) {

  const _a = Math.abs(a);
  const _b = Math.abs(b);

  if (a === 0) {
    return Math.log(_b);
  }

  if (b === 0) {
    return Math.log(_a);
  }

  if (_a < 3000 && _b < 3000) {
    return 0.5 * Math.log(a * a + b * b);
  }

  a = a / 2;
  b = b / 2;

  return 0.5 * Math.log(a * a + b * b) + Math.LN2;
}

/*
 * Temporary parsing object to avoid re-allocations
 *
 */
const P = Object.create(Quaternion.prototype);

function parse(dest, w, x, y, z) {

  // Most common internal use case with 4 params
  if (z !== undefined) {
    dest['w'] = w;
    dest['x'] = x;
    dest['y'] = y;
    dest['z'] = z;
    return;
  }

  if (typeof w === 'object' && y === undefined) {

    // Check for quats, for example when an object gets cloned
    if ('w' in w || 'x' in w || 'y' in w || 'z' in w) {
      dest['w'] = w['w'] || 0;
      dest['x'] = w['x'] || 0;
      dest['y'] = w['y'] || 0;
      dest['z'] = w['z'] || 0;
      return;
    }

    // Check for complex numbers
    if ('re' in w && 'im' in w) {
      dest['w'] = w['re'];
      dest['x'] = w['im'];
      dest['y'] = 0;
      dest['z'] = 0;
      return;
    }

    // Check for array
    if (w.length === 4) {
      dest['w'] = w[0];
      dest['x'] = w[1];
      dest['y'] = w[2];
      dest['z'] = w[3];
      return;
    }

    // Check for augmented vector
    if (w.length === 3) {
      dest['w'] = 0;
      dest['x'] = w[0];
      dest['y'] = w[1];
      dest['z'] = w[2];
      return;
    }

    throw new Error('Invalid object');
  }

  // Parse string values
  if (typeof w === 'string' && y === undefined) {

    const tokens = w.match(/\d+\.?\d*e[+-]?\d+|\d+\.?\d*|\.\d+|./g);
    let plus = 1;
    let minus = 0;

    const iMap = { 'i': 'x', 'j': 'y', 'k': 'z' };

    if (tokens === null) {
      throw new Error('Parse error');
    }

    // Reset the current state
    dest['w'] =
      dest['x'] =
      dest['y'] =
      dest['z'] = 0;

    for (let i = 0; i < tokens.length; i++) {

      let c = tokens[i];
      let d = tokens[i + 1];

      if (c === ' ' || c === '\t' || c === '\n') {
        /* void */
      } else if (c === '+') {
        plus++;
      } else if (c === '-') {
        minus++;
      } else {

        if (plus + minus === 0) {
          throw new Error('Parse error' + c);
        }
        let g = iMap[c];

        // Is the current token an imaginary sign?
        if (g !== undefined) {

          // Is the following token a number?
          if (d !== ' ' && !isNaN(d)) {
            c = d;
            i++;
          } else {
            c = '1';
          }

        } else {

          if (isNaN(c)) {
            throw new Error('Parser error');
          }

          g = iMap[d];

          if (g !== undefined) {
            i++;
          }
        }

        dest[g || 'w'] += parseFloat((minus % 2 ? '-' : '') + c);
        plus = minus = 0;
      }
    }

    // Still something on the stack
    if (plus + minus > 0) {
      throw new Error('Parser error');
    }
    return;
  }

  // If no single variable was given AND it was the constructor, set it to the identity
  if (w === undefined && dest !== P) {
    dest['w'] = 1;
    dest['x'] =
      dest['y'] =
      dest['z'] = 0;
  } else {

    dest['w'] = w || 0;

    // Note: This isn't fromAxis(), it's just syntactic sugar!
    if (x && x.length === 3) {
      dest['x'] = x[0];
      dest['y'] = x[1];
      dest['z'] = x[2];
    } else {
      dest['x'] = x || 0;
      dest['y'] = y || 0;
      dest['z'] = z || 0;
    }
  }
}

function numToStr(n, char, prev) {

  let ret = '';

  if (n !== 0) {

    if (prev !== '') {
      ret += n < 0 ? ' - ' : ' + ';
    } else if (n < 0) {
      ret += '-';
    }

    n = Math.abs(n);

    if (1 !== n || char === '') {
      ret += n;
    }
    ret += char;
  }
  return ret;
}

/**
 * Quaternion constructor
 *
 * @constructor
 * @param {number|Object|string} w real
 * @param {number=} x imag
 * @param {number=} y imag
 * @param {number=} z imag
 * @returns {Quaternion}
 */
function Quaternion(w, x, y, z) {

  if (this instanceof Quaternion) {
    parse(this, w, x, y, z);
  } else {
    const t = Object.create(Quaternion.prototype);
    parse(t, w, x, y, z);
    return t;
  }
}

Quaternion.prototype = {
  'w': 1,
  'x': 0,
  'y': 0,
  'z': 0,
  /**
   * Adds two quaternions Q1 and Q2
   *
   * @param {number|Object|string} w real
   * @param {number=} x imag
   * @param {number=} y imag
   * @param {number=} z imag
   * @returns {Quaternion}
   */
  'add': function (w, x, y, z) {

    parse(P, w, x, y, z);

    // Q1 + Q2 := [w1, v1] + [w2, v2] = [w1 + w2, v1 + v2]

    return newQuaternion(
      this['w'] + P['w'],
      this['x'] + P['x'],
      this['y'] + P['y'],
      this['z'] + P['z']);
  },
  /**
   * Subtracts a quaternions Q2 from Q1
   *
   * @param {number|Object|string} w real
   * @param {number=} x imag
   * @param {number=} y imag
   * @param {number=} z imag
   * @returns {Quaternion}
   */
  'sub': function (w, x, y, z) {

    parse(P, w, x, y, z);

    // Q1 - Q2 := Q1 + (-Q2)
    //          = [w1, v1] - [w2, v2] = [w1 - w2, v1 - v2]

    return newQuaternion(
      this['w'] - P['w'],
      this['x'] - P['x'],
      this['y'] - P['y'],
      this['z'] - P['z']);
  },
  /**
   * Calculates the additive inverse, or simply it negates the quaternion
   *
   * @returns {Quaternion}
   */
  'neg': function () {

    // -Q := [-w, -v]

    return newQuaternion(-this['w'], -this['x'], -this['y'], -this['z']);
  },
  /**
   * Calculates the length/modulus/magnitude or the norm of a quaternion
   *
   * @returns {number}
   */
  'norm': function () {

    // |Q| := sqrt(|Q|^2)

    // The unit quaternion has |Q| = 1

    const w = this['w'];
    const x = this['x'];
    const y = this['y'];
    const z = this['z'];

    return Math.sqrt(w * w + x * x + y * y + z * z);
  },
  /**
   * Calculates the squared length/modulus/magnitude or the norm of a quaternion
   *
   * @returns {number}
   */
  'normSq': function () {

    // |Q|^2 := [w, v] * [w, -v]
    //        = [w^2 + dot(v, v), -w * v + w * v + cross(v, -v)]
    //        = [w^2 + |v|^2, 0]
    //        = [w^2 + dot(v, v), 0]
    //        = dot(Q, Q)
    //        = Q * Q'

    const w = this['w'];
    const x = this['x'];
    const y = this['y'];
    const z = this['z'];

    return w * w + x * x + y * y + z * z;
  },
  /**
   * Normalizes the quaternion to have |Q| = 1 as long as the norm is not zero
   * Alternative names are the signum, unit or versor
   *
   * @returns {Quaternion}
   */
  'normalize': function () {

    // Q* := Q / |Q|

    // unrolled Q.scale(1 / Q.norm())

    const w = this['w'];
    const x = this['x'];
    const y = this['y'];
    const z = this['z'];

    let norm = Math.sqrt(w * w + x * x + y * y + z * z);

    if (norm < EPSILON) {
      return Quaternion['ZERO'];
    }

    norm = 1 / norm;

    return newQuaternion(w * norm, x * norm, y * norm, z * norm);
  },
  /**
   * Calculates the Hamilton product of two quaternions
   * Leaving out the imaginary part results in just scaling the quat
   *
   * @param {number|Object|string} w real
   * @param {number=} x imag
   * @param {number=} y imag
   * @param {number=} z imag
   * @returns {Quaternion}
   */
  'mul': function (w, x, y, z) {

    parse(P, w, x, y, z);

    // Q1 * Q2 = [w1 * w2 - dot(v1, v2), w1 * v2 + w2 * v1 + cross(v1, v2)]

    // Not commutative because cross(v1, v2) != cross(v2, v1)!

    const w1 = this['w'];
    const x1 = this['x'];
    const y1 = this['y'];
    const z1 = this['z'];

    const w2 = P['w'];
    const x2 = P['x'];
    const y2 = P['y'];
    const z2 = P['z'];

    return newQuaternion(
      w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
      w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
      w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2,
      w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2);
  },
  /**
   * Scales a quaternion by a scalar, faster than using multiplication
   *
   * @param {number} s scaling factor
   * @returns {Quaternion}
   */
  'scale': function (s) {

    return newQuaternion(
      this['w'] * s,
      this['x'] * s,
      this['y'] * s,
      this['z'] * s);
  },
  /**
   * Calculates the dot product of two quaternions
   *
   * @param {number|Object|string} w real
   * @param {number=} x imag
   * @param {number=} y imag
   * @param {number=} z imag
   * @returns {number}
   */
  'dot': function (w, x, y, z) {

    parse(P, w, x, y, z);

    // dot(Q1, Q2) := w1 * w2 + dot(v1, v2)

    return this['w'] * P['w'] + this['x'] * P['x'] + this['y'] * P['y'] + this['z'] * P['z'];
  },
  /**
   * Calculates the inverse of a quat for non-normalized quats such that
   * Q^-1 * Q = 1 and Q * Q^-1 = 1
   *
   * @returns {Quaternion}
   */
  'inverse': function () {

    // Q^-1 := Q' / |Q|^2
    //       = [w / (w^2 + |v|^2), -v / (w^2 + |v|^2)]

    // Proof:
    // Q * Q^-1 = [w, v] * [w / (w^2 + |v|^2), -v / (w^2 + |v|^2)]
    //          = [1, 0]
    // Q^-1 * Q = [w / (w^2 + |v|^2), -v / (w^2 + |v|^2)] * [w, v]
    //          = [1, 0].

    const w = this['w'];
    const x = this['x'];
    const y = this['y'];
    const z = this['z'];

    let normSq = w * w + x * x + y * y + z * z;

    if (normSq === 0) {
      return Quaternion['ZERO']; // TODO: Is the result zero or one when the norm is zero?
    }

    normSq = 1 / normSq;

    return newQuaternion(w * normSq, -x * normSq, -y * normSq, -z * normSq);
  },
  /**
   * Multiplies a quaternion with the inverse of a second quaternion
   *
   * @param {number|Object|string} w real
   * @param {number=} x imag
   * @param {number=} y imag
   * @param {number=} z imag
   * @returns {Quaternion}
   */
  'div': function (w, x, y, z) {

    parse(P, w, x, y, z);

    // Q1 / Q2 := Q1 * Q2^-1

    const w1 = this['w'];
    const x1 = this['x'];
    const y1 = this['y'];
    const z1 = this['z'];

    const w2 = P['w'];
    const x2 = P['x'];
    const y2 = P['y'];
    const z2 = P['z'];

    let normSq = w2 * w2 + x2 * x2 + y2 * y2 + z2 * z2;

    if (normSq === 0) {
      return Quaternion['ZERO']; // TODO: Is the result zero or one when the norm is zero?
    }

    normSq = 1 / normSq;

    return newQuaternion(
      (w1 * w2 + x1 * x2 + y1 * y2 + z1 * z2) * normSq,
      (x1 * w2 - w1 * x2 - y1 * z2 + z1 * y2) * normSq,
      (y1 * w2 - w1 * y2 - z1 * x2 + x1 * z2) * normSq,
      (z1 * w2 - w1 * z2 - x1 * y2 + y1 * x2) * normSq);
  },
  /**
   * Calculates the conjugate of a quaternion
   *
   * @returns {Quaternion}
   */
  'conjugate': function () {

    // Q' = [s, -v]

    // If the quaternion is normalized,
    // the conjugate is the inverse of the quaternion - but faster
    // Q' * Q = Q * Q' = 1

    // Additionally, the conjugate of a unit quaternion is a rotation with the same
    // angle but the opposite axis.

    // Moreover the following property holds:
    // (Q1 * Q2)' = Q2' * Q1'

    return newQuaternion(this['w'], -this['x'], -this['y'], -this['z']);
  },
  /**
   * Calculates the natural exponentiation of the quaternion
   *
   * @returns {Quaternion}
   */
  'exp': function () {

    const w = this['w'];
    const x = this['x'];
    const y = this['y'];
    const z = this['z'];

    const vNorm = Math.sqrt(x * x + y * y + z * z);
    const wExp = Math.exp(w);
    const scale = wExp * Math.sin(vNorm) / vNorm;

    if (vNorm === 0) {
      //return newQuaternion(wExp * Math.cos(vNorm), 0, 0, 0);
      return newQuaternion(wExp, 0, 0, 0);
    }

    return newQuaternion(
      wExp * Math.cos(vNorm),
      x * scale,
      y * scale,
      z * scale);
  },
  /**
   * Calculates the natural logarithm of the quaternion
   *
   * @returns {Quaternion}
   */
  'log': function () {

    const w = this['w'];
    const x = this['x'];
    const y = this['y'];
    const z = this['z'];

    if (y === 0 && z === 0) {
      return newQuaternion(
        logHypot(w, x),
        Math.atan2(x, w), 0, 0);
    }

    const qNorm2 = x * x + y * y + z * z + w * w;
    const vNorm = Math.sqrt(x * x + y * y + z * z);

    const scale = Math.atan2(vNorm, w) / vNorm; // Alternative: acos(w / qNorm) / vNorm

    return newQuaternion(
      Math.log(qNorm2) * 0.5,
      x * scale,
      y * scale,
      z * scale);
  },
  /**
   * Calculates the power of a quaternion raised to a real number or another quaternion
   *
   * @param {number|Object|string} w real
   * @param {number=} x imag
   * @param {number=} y imag
   * @param {number=} z imag
   * @returns {Quaternion}
   */
  'pow': function (w, x, y, z) {

    parse(P, w, x, y, z);

    if (P['y'] === 0 && P['z'] === 0) {

      if (P['w'] === 1 && P['x'] === 0) {
        return this;
      }

      if (P['w'] === 0 && P['x'] === 0) {
        return Quaternion['ONE'];
      }

      // Check if we can operate in C
      // Borrowed from complex.js
      if (this['y'] === 0 && this['z'] === 0) {

        let a = this['w'];
        let b = this['x'];

        if (a === 0 && b === 0) {
          return Quaternion['ZERO'];
        }

        let arg = Math.atan2(b, a);
        let loh = logHypot(a, b);

        if (P['x'] === 0) {

          if (b === 0 && a >= 0) {

            return newQuaternion(Math.pow(a, P['w']), 0, 0, 0);

          } else if (a === 0) {

            switch (P['w'] % 4) {
              case 0:
                return newQuaternion(Math.pow(b, P['w']), 0, 0, 0);
              case 1:
                return newQuaternion(0, Math.pow(b, P['w']), 0, 0);
              case 2:
                return newQuaternion(-Math.pow(b, P['w']), 0, 0, 0);
              case 3:
                return newQuaternion(0, -Math.pow(b, P['w']), 0, 0);
            }
          }
        }

        a = Math.exp(P['w'] * loh - P['x'] * arg);
        b = P['x'] * loh + P['w'] * arg;
        return newQuaternion(
          a * Math.cos(b),
          a * Math.sin(b), 0, 0);
      }
    }

    // Normal quaternion behavior
    // q^p = e^ln(q^p) = e^(ln(q)*p)
    return this['log']()['mul'](P)['exp']();
  },
  /**
   * Checks if two quats are the same
   *
   * @param {number|Object|string} w real
   * @param {number=} x imag
   * @param {number=} y imag
   * @param {number=} z imag
   * @returns {boolean}
   */
  'equals': function (w, x, y, z) {

    parse(P, w, x, y, z);

    const eps = EPSILON;

    // maybe check for NaN's here?
    return Math.abs(P['w'] - this['w']) < eps
      && Math.abs(P['x'] - this['x']) < eps
      && Math.abs(P['y'] - this['y']) < eps
      && Math.abs(P['z'] - this['z']) < eps;
  },
  /**
   * Checks if all parts of a quaternion are finite
   *
   * @returns {boolean}
   */
  'isFinite': function () {

    return isFinite(this['w']) && isFinite(this['x']) && isFinite(this['y']) && isFinite(this['z']);
  },
  /**
   * Checks if any of the parts of the quaternion is not a number
   *
   * @returns {boolean}
   */
  'isNaN': function () {

    return isNaN(this['w']) || isNaN(this['x']) || isNaN(this['y']) || isNaN(this['z']);
  },
  /**
   * Gets the Quaternion as a well formatted string
   *
   * @returns {string}
   */
  'toString': function () {

    const w = this['w'];
    const x = this['x'];
    const y = this['y'];
    const z = this['z'];
    let ret = '';

    if (isNaN(w) || isNaN(x) || isNaN(y) || isNaN(z)) {
      return 'NaN';
    }

    // Alternative design?
    // '(%f, [%f %f %f])'

    ret = numToStr(w, '', ret);
    ret += numToStr(x, 'i', ret);
    ret += numToStr(y, 'j', ret);
    ret += numToStr(z, 'k', ret);

    if ('' === ret)
      return '0';

    return ret;
  },
  /**
   * Returns the real part of the quaternion
   *
   * @returns {number}
   */
  'real': function () {

    return this['w'];
  },
  /**
   * Returns the imaginary part of the quaternion as a 3D vector / array
   *
   * @returns {Array}
   */
  'imag': function () {

    return [this['x'], this['y'], this['z']];
  },
  /**
   * Gets the actual quaternion as a 4D vector / array
   *
   * @returns {Array}
   */
  'toVector': function () {

    return [this['w'], this['x'], this['y'], this['z']];
  },
  /**
   * Calculates the 3x3 rotation matrix for the current quat
   *
   * @param {boolean=} twoD
   * @returns {Array}
   */
  'toMatrix': function (twoD) {

    const w = this['w'];
    const x = this['x'];
    const y = this['y'];
    const z = this['z'];

    const wx = w * x, wy = w * y, wz = w * z;
    const xx = x * x, xy = x * y, xz = x * z;
    const yy = y * y, yz = y * z, zz = z * z;

    if (twoD) {
      return [
        [1 - 2 * (yy + zz), 2 * (xy - wz), 2 * (xz + wy)],
        [2 * (xy + wz), 1 - 2 * (xx + zz), 2 * (yz - wx)],
        [2 * (xz - wy), 2 * (yz + wx), 1 - 2 * (xx + yy)]];
    }

    return [
      1 - 2 * (yy + zz), 2 * (xy - wz), 2 * (xz + wy),
      2 * (xy + wz), 1 - 2 * (xx + zz), 2 * (yz - wx),
      2 * (xz - wy), 2 * (yz + wx), 1 - 2 * (xx + yy)];
  },
  /**
   * Calculates the homogeneous 4x4 rotation matrix for the current quat
   *
   * @param {boolean=} twoD
   * @returns {Array}
   */
  'toMatrix4': function (twoD) {

    const w = this['w'];
    const x = this['x'];
    const y = this['y'];
    const z = this['z'];

    const wx = w * x, wy = w * y, wz = w * z;
    const xx = x * x, xy = x * y, xz = x * z;
    const yy = y * y, yz = y * z, zz = z * z;

    if (twoD) {
      return [
        [1 - 2 * (yy + zz), 2 * (xy - wz), 2 * (xz + wy), 0],
        [2 * (xy + wz), 1 - 2 * (xx + zz), 2 * (yz - wx), 0],
        [2 * (xz - wy), 2 * (yz + wx), 1 - 2 * (xx + yy), 0],
        [0, 0, 0, 1]];
    }

    return [
      1 - 2 * (yy + zz), 2 * (xy - wz), 2 * (xz + wy), 0,
      2 * (xy + wz), 1 - 2 * (xx + zz), 2 * (yz - wx), 0,
      2 * (xz - wy), 2 * (yz + wx), 1 - 2 * (xx + yy), 0,
      0, 0, 0, 1];
  },
  /**
   * Determines the homogeneous rotation matrix string used for CSS 3D transforms
   *
   * @returns {string}
   */
  'toCSSTransform': function () {

    const w = this['w'];

    let angle = 2 * Math.acos(w);
    let sin2 = 1 - w * w;

    if (sin2 < EPSILON) {
      angle = 0;
      sin2 = 1;
    } else {
      sin2 = 1 / Math.sqrt(sin2); // Re-use variable sin^2 for 1 / sin
    }
    return "rotate3d(" + this['x'] * sin2 + "," + this['y'] * sin2 + "," + this['z'] * sin2 + "," + angle + "rad)";
  },
  /**
   * Calculates the axis and angle representation of the quaternion
   *
   * @returns {Array}
   */
  'toAxisAngle': function () {

    const w = this['w'];
    const sin2 = 1 - w * w; // sin(angle / 2) = sin(acos(w)) = sqrt(1 - w^2) = |v|, since 1 = dot(Q) <=> dot(v) = 1 - w^2

    if (sin2 < EPSILON) { // Alternatively |v| == 0
      // If the sine is close to 0, we're close to the unit quaternion and the direction does not matter
      return [[this['x'], this['y'], this['z']], 0]; // or [[1, 0, 0], 0] ?  or [[0, 0, 0], 0] ?
    }

    const isin = 1 / Math.sqrt(sin2);
    const angle = 2 * Math.acos(w); // Alternatively: 2 * atan2(|v|, w)
    return [[this['x'] * isin, this['y'] * isin, this['z'] * isin], angle];
  },
  /**
   * Calculates the Euler angles represented by the current quat (multiplication order from right to left)
   * 
   * @param {string=} order Axis order (Tait Bryan)
   * @returns {Object}
   */
  'toEuler': function (order) {

    const w = this['w'];
    const x = this['x'];
    const y = this['y'];
    const z = this['z'];

    const wx = w * x, wy = w * y, wz = w * z;
    const xx = x * x, xy = x * y, xz = x * z;
    const yy = y * y, yz = y * z, zz = z * z;

    function asin(t) {
      return t >= 1 ? Math.PI / 2 : (t <= -1 ? -Math.PI / 2 : Math.asin(t));
    }

    if (order === undefined || order === 'ZXY') {
      return [
        -Math.atan2(2 * (xy - wz), 1 - 2 * (xx + zz)),
        asin(2 * (yz + wx)),
        -Math.atan2(2 * (xz - wy), 1 - 2 * (xx + yy)),
      ];
    }

    if (order === 'XYZ' || order === 'RPY') {
      return [
        -Math.atan2(2 * (yz - wx), 1 - 2 * (xx + yy)),
        asin(2 * (xz + wy)),
        -Math.atan2(2 * (xy - wz), 1 - 2 * (yy + zz)),
      ];
    }

    if (order === 'YXZ') {
      return [
        Math.atan2(2 * (xz + wy), 1 - 2 * (xx + yy)),
        -asin(2 * (yz - wx)),
        Math.atan2(2 * (xy + wz), 1 - 2 * (xx + zz)),
      ];
    }

    if (order === 'ZYX' || order === 'YPR') {  // roll around X, pitch around Y, yaw around Z
      /*
      if (2 * (xz - wy) > .999) {
        return [
          2 * Math.atan2(x, w),
          -Math.PI / 2,
          0
        ];
      }

      if (2 * (xz - wy) < -.999) {
        return [
          -2 * Math.atan2(x, w),
          Math.PI / 2,
          0
        ];
      }
      */
      return [
        Math.atan2(2 * (xy + wz), 1 - 2 * (yy + zz)), // Heading / Yaw
        -asin(2 * (xz - wy)), // Attitude / Pitch
        Math.atan2(2 * (yz + wx), 1 - 2 * (xx + yy)), // Bank / Roll
      ];
    }

    if (order === 'YZX') {
      /*
      if (2 * (xy + wz) > .999) { // North pole
        return [
          2 * Math.atan2(x, w),
          Math.PI / 2,
          0
        ];
      }

      if (2 * (xy + wz) < -.999) { // South pole
        return [
          -2 * Math.atan2(x, w),
          -Math.PI / 2,
          0
        ];
      }
      */
      return [
        -Math.atan2(2 * (xz - wy), 1 - 2 * (yy + zz)), // Heading
        asin(2 * (xy + wz)), // Attitude
        -Math.atan2(2 * (yz - wx), 1 - 2 * (xx + zz)), // Bank
      ];
    }

    if (order === 'XZY') {
      return [
        Math.atan2(2 * (yz + wx), 1 - 2 * (xx + zz)),
        -asin(2 * (xy - wz)),
        Math.atan2(2 * (xz + wy), 1 - 2 * (yy + zz)),
      ];
    }
    return null;
  },
  /**
   * Clones the actual object
   *
   * @returns {Quaternion}
   */
  'clone': function () {

    return newQuaternion(this['w'], this['x'], this['y'], this['z']);
  },
  /**
   * Rotates a vector according to the current quaternion, assumes |q|=1
   * @link https://raw.org/proof/vector-rotation-using-quaternions/
   *
   * @param {Array} v The vector to be rotated
   * @returns {Array}
   */
  'rotateVector': function (v) {

    const qw = this['w'];
    const qx = this['x'];
    const qy = this['y'];
    const qz = this['z'];

    const vx = v[0];
    const vy = v[1];
    const vz = v[2];

    // t = q x v
    let tx = qy * vz - qz * vy;
    let ty = qz * vx - qx * vz;
    let tz = qx * vy - qy * vx;

    // t = 2t
    tx = tx + tx;
    ty = ty + ty;
    tz = tz + tz;

    // v + w t + q x t
    return [
      vx + qw * tx + qy * tz - qz * ty,
      vy + qw * ty + qz * tx - qx * tz,
      vz + qw * tz + qx * ty - qy * tx];
  },

  /**
   * Gets a function to spherically interpolate between two quaternions
   * 
   * @returns Function
   */
  'slerp': function (w, x, y, z) {

    parse(P, w, x, y, z);

    // slerp(Q1, Q2, t) := Q1(Q1^-1 Q2)^t

    let w1 = this['w'];
    let x1 = this['x'];
    let y1 = this['y'];
    let z1 = this['z'];

    let w2 = P['w'];
    let x2 = P['x'];
    let y2 = P['y'];
    let z2 = P['z'];

    let cosTheta0 = w1 * w2 + x1 * x2 + y1 * y2 + z1 * z2;

    if (cosTheta0 < 0) {
      w1 = -w1;
      x1 = -x1;
      y1 = -y1;
      z1 = -z1;
      cosTheta0 = -cosTheta0;
    }

    if (cosTheta0 >= 1 - EPSILON) {
      return function (pct) {
        return newNormalized(
          w1 + pct * (w2 - w1),
          x1 + pct * (x2 - x1),
          y1 + pct * (y2 - y1),
          z1 + pct * (z2 - z1));
      };
    }

    let Theta0 = Math.acos(cosTheta0);
    let sinTheta0 = Math.sin(Theta0);

    return function (pct) {

      let Theta = Theta0 * pct;
      let sinTheta = Math.sin(Theta);
      let cosTheta = Math.cos(Theta);

      let s0 = cosTheta - cosTheta0 * sinTheta / sinTheta0;
      let s1 = sinTheta / sinTheta0;

      return newQuaternion(
        s0 * w1 + s1 * w2,
        s0 * x1 + s1 * x2,
        s0 * y1 + s1 * y2,
        s0 * z1 + s1 * z2);
    };
  }
};

Quaternion['ZERO'] = newQuaternion(0, 0, 0, 0); // This is the additive identity Quaternion
Quaternion['ONE'] = newQuaternion(1, 0, 0, 0); // This is the multiplicative identity Quaternion
Quaternion['I'] = newQuaternion(0, 1, 0, 0);
Quaternion['J'] = newQuaternion(0, 0, 1, 0);
Quaternion['K'] = newQuaternion(0, 0, 0, 1);

/**
 * @const
 */
const EPSILON = 1e-16;

/**
 * Creates quaternion by a rotation given as axis-angle orientation
 *
 * @param {Array} axis The axis around which to rotate
 * @param {number} angle The angle in radians
 * @returns {Quaternion}
 */
Quaternion['fromAxisAngle'] = function (axis, angle) {

  // Q = [cos(angle / 2), v * sin(angle / 2)]

  const a = axis[0];
  const b = axis[1];
  const c = axis[2];

  const halfAngle = angle * 0.5;

  const sin_2 = Math.sin(halfAngle);
  const cos_2 = Math.cos(halfAngle);

  const sin_norm = sin_2 / Math.sqrt(a * a + b * b + c * c);

  return newQuaternion(cos_2, a * sin_norm, b * sin_norm, c * sin_norm);
};

/**
 * Calculates the quaternion to rotate vector u onto vector v
 * @link https://raw.org/proof/quaternion-from-two-vectors/
 *
 * @param {Array} u
 * @param {Array} v
 */
Quaternion['fromVectors'] = function (u, v) {

  let ux = u[0];
  let uy = u[1];
  let uz = u[2];

  let vx = v[0];
  let vy = v[1];
  let vz = v[2];

  const uLen = Math.sqrt(ux * ux + uy * uy + uz * uz);
  const vLen = Math.sqrt(vx * vx + vy * vy + vz * vz);

  // Normalize u and v
  if (uLen > 0) ux /= uLen, uy /= uLen, uz /= uLen;
  if (vLen > 0) vx /= vLen, vy /= vLen, vz /= vLen;

  // Calculate dot product of normalized u and v
  const dot = ux * vx + uy * vy + uz * vz;

  // Parallel when dot > 0.999999
  if (dot >= 1 - EPSILON) {
    return Quaternion['ONE'];
  }

  // Anti-Parallel (close to PI) when dot < -0.999999
  if (1 + dot <= EPSILON) {

    // Rotate 180° around any orthogonal vector
    // axis = len(cross([1, 0, 0], u)) == 0 ? cross([0, 1, 0], u) : cross([1, 0, 0], u) and therefore
    //    return Quaternion['fromAxisAngle'](Math.abs(ux) > Math.abs(uz) ? [-uy, ux, 0] : [0, -uz, uy], Math.PI)
    // or return Quaternion['fromAxisAngle'](Math.abs(ux) > Math.abs(uz) ? [ uy,-ux, 0] : [0,  uz,-uy], Math.PI)
    // or ...

    // Since fromAxisAngle(axis, PI) == Quaternion(0, axis).normalize(),
    if (Math.abs(ux) > Math.abs(uz)) {
      return newNormalized(0, -uy, ux, 0);
    } else {
      return newNormalized(0, 0, -uz, uy);
    }
  }

  // w = cross(u, v)
  const wx = uy * vz - uz * vy;
  const wy = uz * vx - ux * vz;
  const wz = ux * vy - uy * vx;

  // |Q| = sqrt((1.0 + dot) * 2.0)
  return newNormalized(1 + dot, wx, wy, wz);
};

/**
 * Gets a spherical random number
 * @link http://planning.cs.uiuc.edu/node198.html
 */
Quaternion['random'] = function () {

  const u1 = Math.random();
  const u2 = Math.random();
  const u3 = Math.random();

  const s = Math.sqrt(1 - u1);
  const t = Math.sqrt(u1);

  return newQuaternion(
    t * Math.cos(2 * Math.PI * u3),
    s * Math.sin(2 * Math.PI * u2),
    s * Math.cos(2 * Math.PI * u2),
    t * Math.sin(2 * Math.PI * u3)
  );
};

/**
 * Creates a quaternion by a rotation given by Euler angles (logical application order from left to right)
 *
 * @param {number} ψ First angle
 * @param {number} θ Second angle
 * @param {number} φ Third angle
 * @param {string=} order Axis order (Tait Bryan)
 * @returns {Quaternion}
 */
Quaternion['fromEulerLogical'] = function (ψ, θ, φ, order) {

  return Quaternion['fromEuler'](φ, θ, ψ, order !== undefined ? order[2] + order[1] + order[0] : order);
};

/**
 * Creates a quaternion by a rotation given by Euler angles (multiplication order from right to left)
 *
 * @param {number} φ First angle
 * @param {number} θ Second angle
 * @param {number} ψ Third angle
 * @param {string=} order Axis order (Tait Bryan)
 * @returns {Quaternion}
 */
Quaternion['fromEuler'] = function (φ, θ, ψ, order) {

  const _x = φ * 0.5;
  const _y = θ * 0.5;
  const _z = ψ * 0.5;

  const cX = Math.cos(_x);
  const cY = Math.cos(_y);
  const cZ = Math.cos(_z);

  const sX = Math.sin(_x);
  const sY = Math.sin(_y);
  const sZ = Math.sin(_z);

  if (order === undefined || order === 'ZXY') {
    // axisAngle([0, 0, 1], φ) * axisAngle([1, 0, 0], θ) * axisAngle([0, 1, 0], ψ)
    return newQuaternion(
      cX * cY * cZ - sX * sY * sZ,
      sY * cX * cZ - sX * sZ * cY,
      sX * sY * cZ + sZ * cX * cY,
      sX * cY * cZ + sY * sZ * cX);
  }

  if (order === 'XYZ' || order === 'RPY') { // roll around X, pitch around Y, yaw around Z
    // axisAngle([1, 0, 0], φ) * axisAngle([0, 1, 0], θ) * axisAngle([0, 0, 1], ψ)
    return newQuaternion(
      cX * cY * cZ - sX * sY * sZ,
      sX * cY * cZ + sY * sZ * cX,
      sY * cX * cZ - sX * sZ * cY,
      sX * sY * cZ + sZ * cX * cY);
  }

  if (order === 'YXZ') { // deviceorientation
    // axisAngle([0, 1, 0], φ) * axisAngle([1, 0, 0], θ) * axisAngle([0, 0, 1], ψ)
    return newQuaternion(
      sX * sY * sZ + cX * cY * cZ,
      sX * sZ * cY + sY * cX * cZ,
      sX * cY * cZ - sY * sZ * cX,
      sZ * cX * cY - sX * sY * cZ);
  }

  if (order === 'ZYX' || order === 'YPR') {
    // axisAngle([0, 0, 1], φ) * axisAngle([0, 1, 0], θ) * axisAngle([1, 0, 0], ψ)
    return newQuaternion(
      sX * sY * sZ + cX * cY * cZ,
      sZ * cX * cY - sX * sY * cZ,
      sX * sZ * cY + sY * cX * cZ,
      sX * cY * cZ - sY * sZ * cX);
  }

  if (order === 'YZX') {
    // axisAngle([0, 1, 0], φ) * axisAngle([0, 0, 1], θ) * axisAngle([1, 0, 0], ψ)
    return newQuaternion(
      cX * cY * cZ - sX * sY * sZ,
      sX * sY * cZ + sZ * cX * cY,
      sX * cY * cZ + sY * sZ * cX,
      sY * cX * cZ - sX * sZ * cY);
  }

  if (order === 'XZY') {
    // axisAngle([1, 0, 0], φ) * axisAngle([0, 0, 1], θ) * axisAngle([0, 1, 0], ψ)
    return newQuaternion(
      sX * sY * sZ + cX * cY * cZ,
      sX * cY * cZ - sY * sZ * cX,
      sZ * cX * cY - sX * sY * cZ,
      sX * sZ * cY + sY * cX * cZ);
  }

  if (order === 'ZYZ') {
    // axisAngle([0, 0, 1], φ) * axisAngle([0, 1, 0], θ) * axisAngle([0, 0, 1], ψ)
    return newQuaternion(
      cX * cY * cZ - sX * sZ * cY,
      sY * sZ * cX - sX * sY * cZ,
      sX * sY * sZ + sY * cX * cZ,
      sX * cY * cZ + sZ * cX * cY);
  }

  if (order === 'ZXZ') {
    // axisAngle([0, 0, 1], φ) * axisAngle([1, 0, 0], θ) * axisAngle([0, 0, 1], ψ)
    return newQuaternion(
      cX * cY * cZ - sX * sZ * cY,
      sX * sY * sZ + sY * cX * cZ,
      sX * sY * cZ - sY * sZ * cX,
      sX * cY * cZ + sZ * cX * cY);
  }

  if (order === 'YXY') {
    // axisAngle([0, 1, 0], φ) * axisAngle([1, 0, 0], θ) * axisAngle([0, 1, 0], ψ)
    return newQuaternion(
      cX * cY * cZ - sX * sZ * cY,
      sX * sY * sZ + sY * cX * cZ,
      sX * cY * cZ + sZ * cX * cY,
      sY * sZ * cX - sX * sY * cZ);
  }

  if (order === 'YZY') {
    // axisAngle([0, 1, 0], φ) * axisAngle([0, 0, 1], θ) * axisAngle([0, 1, 0], ψ)
    return newQuaternion(
      cX * cY * cZ - sX * sZ * cY,
      sX * sY * cZ - sY * sZ * cX,
      sX * cY * cZ + sZ * cX * cY,
      sX * sY * sZ + sY * cX * cZ);
  }

  if (order === 'XYX') {
    // axisAngle([1, 0, 0], φ) * axisAngle([0, 1, 0], θ) * axisAngle([1, 0, 0], ψ)
    return newQuaternion(
      cX * cY * cZ - sX * sZ * cY,
      sX * cY * cZ + sZ * cX * cY,
      sX * sY * sZ + sY * cX * cZ,
      sX * sY * cZ - sY * sZ * cX);
  }

  if (order === 'XZX') {
    // axisAngle([1, 0, 0], φ) * axisAngle([0, 0, 1], θ) * axisAngle([1, 0, 0], ψ)
    return newQuaternion(
      cX * cY * cZ - sX * sZ * cY,
      sX * cY * cZ + sZ * cX * cY,
      sY * sZ * cX - sX * sY * cZ,
      sX * sY * sZ + sY * cX * cZ);
  }
  return null;
};

/**
 * Creates a quaternion by a rotation matrix
 *
 * @param {Array} matrix
 * @returns {Quaternion}
 */
Quaternion['fromMatrix'] = function (matrix) {

  let m00, m01, m02, m10, m11, m12, m20, m21, m22;

  if (matrix.length === 9) {
    m00 = matrix[0];
    m01 = matrix[1];
    m02 = matrix[2];

    m10 = matrix[3];
    m11 = matrix[4];
    m12 = matrix[5];

    m20 = matrix[6];
    m21 = matrix[7];
    m22 = matrix[8];

  } else {
    m00 = matrix[0][0];
    m01 = matrix[0][1];
    m02 = matrix[0][2];

    m10 = matrix[1][0];
    m11 = matrix[1][1];
    m12 = matrix[1][2];

    m20 = matrix[2][0];
    m21 = matrix[2][1];
    m22 = matrix[2][2];
  }

  const tr = m00 + m11 + m22; // 2 * w = sqrt(1 + tr)

  // Choose the element with the biggest value on the diagonal

  if (tr > 0) { // if trace is positive then "w" is biggest component
    // |Q| = 2 * sqrt(1 + tr) = 4w
    return newNormalized(
      tr + 1.0,
      m21 - m12,
      m02 - m20,
      m10 - m01);
  } else if (m00 > m11 && m00 > m22) {
    // |Q| = 2 * sqrt(1.0 + m00 - m11 - m22) = 4x
    return newNormalized(
      m21 - m12,
      1.0 + m00 - m11 - m22,
      m01 + m10,
      m02 + m20);
  } else if (m11 > m22) {
    // |Q| = 2 * sqrt(1.0 + m11 - m00 - m22) = 4y
    return newNormalized(
      m02 - m20,
      m01 + m10,
      1.0 + m11 - m00 - m22,
      m12 + m21);
  } else {
    // |Q| = 2 * sqrt(1.0 + m22 - m00 - m11) = 4z
    return newNormalized(
      m10 - m01,
      m02 + m20,
      m12 + m21,
      1.0 + m22 - m00 - m11);
  }
};

Object.defineProperty(Quaternion, "__esModule", { 'value': true });
Quaternion['default'] = Quaternion;
Quaternion['Quaternion'] = Quaternion;
module['exports'] = Quaternion;


/***/ }),

/***/ "./src/aliases.ts":
/*!************************!*\
  !*** ./src/aliases.ts ***!
  \************************/
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.aliases = void 0;
exports.aliases = {
    "mPelvis": "hip",
    "mSpine1": "mSpine1",
    "mSpine2": "mSpine2",
    "mTorso": "abdomen",
    "mSpine3": "mSpine3",
    "mSpine4": "mSpine4",
    "mChest": "chest",
    "mNeck": "neck",
    "mHead": "head",
    "mSkull": "figureHair",
    "mEyeRight": "mEyeRight",
    "mEyeLeft": "mEyeLeft",
    "mFaceRoot": "mFaceRoot",
    "mFaceEyeAltRight": "mFaceEyeAltRight",
    "mFaceEyeAltLeft": "mFaceEyeAltLeft",
    "mFaceForeheadLeft": "mFaceForeheadLeft",
    "mFaceForeheadRight": "mFaceForeheadRight",
    "mFaceEyebrowOuterLeft": "mFaceEyebrowOuterLeft",
    "mFaceEyebrowCenterLeft": "mFaceEyebrowCenterLeft",
    "mFaceEyebrowInnerLeft": "mFaceEyebrowInnerLeft",
    "mFaceEyebrowOuterRight": "mFaceEyebrowOuterRight",
    "mFaceEyebrowCenterRight": "mFaceEyebrowCenterRight",
    "mFaceEyebrowInnerRight": "mFaceEyebrowInnerRight",
    "mFaceEyeLidUpperLeft": "mFaceEyeLidUpperLeft",
    "mFaceEyeLidLowerLeft": "mFaceEyeLidLowerLeft",
    "mFaceEyeLidUpperRight": "mFaceEyeLidUpperRight",
    "mFaceEyeLidLowerRight": "mFaceEyeLidLowerRight",
    "mFaceEar1Left": "mFaceEar1Left",
    "mFaceEar2Left": "mFaceEar2Left",
    "mFaceEar1Right": "mFaceEar1Right",
    "mFaceEar2Right": "mFaceEar2Right",
    "mFaceNoseLeft": "mFaceNoseLeft",
    "mFaceNoseCenter": "mFaceNoseCenter",
    "mFaceNoseRight": "mFaceNoseRight",
    "mFaceCheekLowerLeft": "mFaceCheekLowerLeft",
    "mFaceCheekUpperLeft": "mFaceCheekUpperLeft",
    "mFaceCheekLowerRight": "mFaceCheekLowerRight",
    "mFaceCheekUpperRight": "mFaceCheekUpperRight",
    "mFaceJaw": "mFaceJaw",
    "mFaceChin": "mFaceChin",
    "mFaceTeethLower": "mFaceTeethLower",
    "mFaceLipLowerLeft": "mFaceLipLowerLeft",
    "mFaceLipLowerRight": "mFaceLipLowerRight",
    "mFaceLipLowerCenter": "mFaceLipLowerCenter",
    "mFaceTongueBase": "mFaceTongueBase",
    "mFaceTongueTip": "mFaceTongueTip",
    "mFaceJawShaper": "mFaceJawShaper",
    "mFaceForeheadCenter": "mFaceForeheadCenter",
    "mFaceNoseBase": "mFaceNoseBase",
    "mFaceTeethUpper": "mFaceTeethUpper",
    "mFaceLipUpperLeft": "mFaceLipUpperLeft",
    "mFaceLipUpperRight": "mFaceLipUpperRight",
    "mFaceLipCornerLeft": "mFaceLipCornerLeft",
    "mFaceLipCornerRight": "mFaceLipCornerRight",
    "mFaceLipUpperCenter": "mFaceLipUpperCenter",
    "mFaceEyecornerInnerLeft": "mFaceEyecornerInnerLeft",
    "mFaceEyecornerInnerRight": "mFaceEyecornerInnerRight",
    "mFaceNoseBridge": "mFaceNoseBridge",
    "mCollarLeft": "lCollar",
    "mShoulderLeft": "lShldr",
    "mElbowLeft": "lForeArm",
    "mWristLeft": "lHand",
    "mHandMiddle1Left": "mHandMiddle1Left",
    "mHandMiddle2Left": "mHandMiddle2Left",
    "mHandMiddle3Left": "mHandMiddle3Left",
    "mHandIndex1Left": "mHandIndex1Left",
    "mHandIndex2Left": "mHandIndex2Left",
    "mHandIndex3Left": "mHandIndex3Left",
    "mHandRing1Left": "mHandRing1Left",
    "mHandRing2Left": "mHandRing2Left",
    "mHandRing3Left": "mHandRing3Left",
    "mHandPinky1Left": "mHandPinky1Left",
    "mHandPinky2Left": "mHandPinky2Left",
    "mHandPinky3Left": "mHandPinky3Left",
    "mHandThumb1Left": "mHandThumb1Left",
    "mHandThumb2Left": "mHandThumb2Left",
    "mHandThumb3Left": "mHandThumb3Left",
    "mCollarRight": "rCollar",
    "mShoulderRight": "rShldr",
    "mElbowRight": "rForeArm",
    "mWristRight": "rHand",
    "mHandMiddle1Right": "mHandMiddle1Right",
    "mHandMiddle2Right": "mHandMiddle2Right",
    "mHandMiddle3Right": "mHandMiddle3Right",
    "mHandIndex1Right": "mHandIndex1Right",
    "mHandIndex2Right": "mHandIndex2Right",
    "mHandIndex3Right": "mHandIndex3Right",
    "mHandRing1Right": "mHandRing1Right",
    "mHandRing2Right": "mHandRing2Right",
    "mHandRing3Right": "mHandRing3Right",
    "mHandPinky1Right": "mHandPinky1Right",
    "mHandPinky2Right": "mHandPinky2Right",
    "mHandPinky3Right": "mHandPinky3Right",
    "mHandThumb1Right": "mHandThumb1Right",
    "mHandThumb2Right": "mHandThumb2Right",
    "mHandThumb3Right": "mHandThumb3Right",
    "mWingsRoot": "mWingsRoot",
    "mWing1Left": "mWing1Left",
    "mWing2Left": "mWing2Left",
    "mWing3Left": "mWing3Left",
    "mWing4Left": "mWing4Left",
    "mWing4FanLeft": "mWing4FanLeft",
    "mWing1Right": "mWing1Right",
    "mWing2Right": "mWing2Right",
    "mWing3Right": "mWing3Right",
    "mWing4Right": "mWing4Right",
    "mWing4FanRight": "mWing4FanRight",
    "mHipRight": "rThigh",
    "mKneeRight": "rShin",
    "mAnkleRight": "rFoot",
    "mFootRight": "mFootRight",
    "mToeRight": "mToeRight",
    "mHipLeft": "lThigh",
    "mKneeLeft": "lShin",
    "mAnkleLeft": "lFoot",
    "mFootLeft": "mFootLeft",
    "mToeLeft": "mToeLeft",
    "mTail1": "mTail1",
    "mTail2": "mTail2",
    "mTail3": "mTail3",
    "mTail4": "mTail4",
    "mTail5": "mTail5",
    "mTail6": "mTail6",
    "mGroin": "mGroin",
    "mHindLimbsRoot": "mHindLimbsRoot",
    "mHindLimb1Left": "mHindLimb1Left",
    "mHindLimb2Left": "mHindLimb2Left",
    "mHindLimb3Left": "mHindLimb3Left",
    "mHindLimb4Left": "mHindLimb4Left",
    "mHindLimb1Right": "mHindLimb1Right",
    "mHindLimb2Right": "mHindLimb2Right",
    "mHindLimb3Right": "mHindLimb3Right",
    "mHindLimb4Right": "mHindLimb4Right"
};


/***/ }),

/***/ "./src/convert.ts":
/*!************************!*\
  !*** ./src/convert.ts ***!
  \************************/
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.visitNode = visitNode;
exports.serializeBVH = serializeBVH;
exports.toBVH = toBVH;
const utils_1 = __webpack_require__(/*! ./utils */ "./src/utils.ts");
const hierarchy_1 = __webpack_require__(/*! ./hierarchy */ "./src/hierarchy.ts");
const aliases_1 = __webpack_require__(/*! ./aliases */ "./src/aliases.ts");
const quaternion_1 = __webpack_require__(/*! quaternion */ "./node_modules/quaternion/dist/quaternion.js");
function offsetToString(offset, digits) {
    return "OFFSET " + (0, utils_1.floatToString)(offset.x, digits) + " " + (0, utils_1.floatToString)(offset.y, digits) + " " + (0, utils_1.floatToString)(offset.z, digits);
}
function channelsString(node) {
    if (!node.channels) {
        return "";
    }
    if (node.bvhName === "hip") {
        return "CHANNELS 6 Xposition Yposition Zposition Xrotation Yrotation Zrotation";
    }
    return "CHANNELS " + node.channels.length + " " + node.channels.join(" ");
}
function appendNode(joint, tabs) {
    let result = "";
    const boneType = (joint.bvhName === "hip") ? "ROOT" : "JOINT";
    const channels = channelsString(joint);
    const offset = (joint.bvhName === "hip") ? offsetToString(joint.offset, 6) : offsetToString(joint.offset, 4);
    if (joint.bvhName != "end") {
        result += tabs + boneType + " " + joint.bvhName + "\n" + tabs + "{\n";
    }
    else {
        result += tabs + "End Site" + "\n" + tabs + "{\n";
    }
    result += tabs + "\t" + offset + "\n";
    if (joint.bvhName != "end") {
        result += tabs + "\t" + channels + "\n";
    }
    if (joint.children) {
        joint.children.forEach((item) => { result += appendNode(item, tabs + "\t"); });
    }
    result += tabs + "}\n";
    return result;
}
function containsNames(node, bvhNames) {
    if (bvhNames.includes(node.bvhName)) {
        return true;
    }
    if (!node.children) {
        return false;
    }
    return !!node.children.map((item) => containsNames(item, bvhNames)).find((item) => !!item);
}
function collectNodes(node, bvhNames) {
    const result = {};
    if (containsNames(node, bvhNames)) {
        result.bvhName = node.bvhName;
    }
    else {
        result.exclude = true;
        return result;
    }
    if (node.children && !!node.children.map((item) => containsNames(item, bvhNames)).find((item) => !!item)) {
        result.children = node.children.map((item) => collectNodes(item, bvhNames)).filter((item) => !item.exclude);
    }
    else {
        result.children = [];
    }
    if (result.children.length > 0) {
        return result;
    }
    result.children.push({ bvhName: "end" });
    return result;
}
function subTree(joints) {
    const names = joints.map(item => item.joint_name);
    const bvhNames = names.map(item => aliases_1.aliases[item] || item);
    return collectNodes(hierarchy_1.hierarchy, bvhNames);
}
function visitNode(node, visitor, childrenFirst = false) {
    if (node.children && childrenFirst) {
        node.children.toReversed().forEach((item) => visitNode(item, visitor, true));
    }
    visitor(node);
    if (node.children && !childrenFirst) {
        node.children.toReversed().forEach((item) => visitNode(item, visitor, false));
    }
}
function extractFramesLength(animJoints) {
    var _a, _b;
    const joint = animJoints.find((item) => { var _a, _b; return ((_a = item.position_keys) === null || _a === void 0 ? void 0 : _a.length) || ((_b = item.rotation_keys) === null || _b === void 0 ? void 0 : _b.length); });
    return ((_a = joint === null || joint === void 0 ? void 0 : joint.position_keys) === null || _a === void 0 ? void 0 : _a.length) || ((_b = joint === null || joint === void 0 ? void 0 : joint.rotation_keys) === null || _b === void 0 ? void 0 : _b.length);
}
function extractTimes(animJoints) {
    var _a;
    const joint = animJoints.find((item) => { var _a, _b; return ((_a = item.position_keys) === null || _a === void 0 ? void 0 : _a.length) || ((_b = item.rotation_keys) === null || _b === void 0 ? void 0 : _b.length); });
    const timeHolders = ((_a = joint === null || joint === void 0 ? void 0 : joint.position_keys) === null || _a === void 0 ? void 0 : _a.length) ? joint.position_keys : joint === null || joint === void 0 ? void 0 : joint.rotation_keys;
    return (timeHolders || []).map((item) => item.time);
}
function fillChannels(node, joint) {
    var _a;
    if (node.bvhName == "hip") {
        node.channels = ["Xposition", "Yposition", "Zposition", "Xrotation", "Yrotation", "Zrotation"];
        return;
    }
    node.channels = [];
    if ((_a = joint === null || joint === void 0 ? void 0 : joint.position_keys) === null || _a === void 0 ? void 0 : _a.length) {
        node.channels.push("Xposition");
        node.channels.push("Yposition");
        node.channels.push("Zposition");
    }
    node.channels.push("Xrotation");
    node.channels.push("Yrotation");
    node.channels.push("Zrotation");
}
function animPositionToBvh(position) {
    const multiplier = 39.3795;
    return { x: position.y * multiplier, y: position.z * multiplier, z: position.x * multiplier };
}
function fillKeyFrames(data, bvhNode, fps) {
    const animJoints = data.joints;
    const length = extractFramesLength(animJoints);
    const bvhTimes = (0, utils_1.getUniformTimes)(data.duration || 0.5, 1 / fps);
    visitNode(bvhNode, (node) => {
        const joint = animJoints.find((item) => aliases_1.aliases[item.joint_name] === node.bvhName);
        node.offset = { x: 0, y: 0, z: 0 };
        if (node.bvhName != "end") {
            fillChannels(node, joint);
            node.animKeys = {
                positions: (joint === null || joint === void 0 ? void 0 : joint.position_keys) || [],
                rotations: (joint === null || joint === void 0 ? void 0 : joint.rotation_keys) || []
            };
            const positionTimes = (0, utils_1.clipTimesToClosestBVHTime)(node.animKeys.positions.map((item) => item.time), bvhTimes);
            const rotationTimes = (0, utils_1.clipTimesToClosestBVHTime)(node.animKeys.rotations.map((item) => item.time), bvhTimes);
            const positions = (0, utils_1.lerpValues)(node.animKeys.positions.map((item) => animPositionToBvh(item)), positionTimes, bvhTimes, { x: 0, y: 0, z: 0 }, utils_1.lerpVector);
            const rotations = (0, utils_1.lerpValues)(node.animKeys.rotations.map((item) => (0, utils_1.toQuaternion)(item)), rotationTimes, bvhTimes, new quaternion_1.Quaternion(), utils_1.lerpQuaternion).map(item => (0, utils_1.quaternionToEulers)(item));
            node.bvhFrames = [];
            bvhTimes.forEach((item, i) => node.bvhFrames.push({
                position: positions[i],
                rotation: rotations[i]
            }));
        }
    });
    bvhNode.bvhTimes = bvhTimes;
}
function getValue(bvhNode, channel, frameNum) {
    const frame = bvhNode.bvhFrames[frameNum];
    const key = channel.toLowerCase()[0];
    const data = (channel.includes("pos")) ? frame.position : frame.rotation;
    const value = data[key];
    return (Math.abs(value) > 0.00000001) ? value : 0;
}
function getValues(bvhNode, frameNum) {
    return bvhNode.channels.map((item) => getValue(bvhNode, item, frameNum));
}
function getFrameValues(bvhNode, frameNum) {
    const result = [];
    visitNode(bvhNode, (node) => {
        if (!node.channels) {
            return;
        }
        result.unshift(...getValues(node, frameNum));
    }, true);
    return result;
}
function getFrameRow(bvhNode, frameNum) {
    const values = getFrameValues(bvhNode, frameNum);
    return values.map(item => (0, utils_1.floatToString)(item, 4)).join(" ") + " \n";
}
function serializeBVH(bvhNode) {
    let result = "HIERARCHY\n";
    result += appendNode(bvhNode, "");
    result += "MOTION\n";
    result += "Frames: " + bvhNode.bvhTimes.length + "\n";
    result += "Frame Time: " + (0, utils_1.floatToString)(bvhNode.bvhTimes[1], 6) + "\n";
    for (let i = 0; i < bvhNode.bvhTimes.length; i++) {
        result += getFrameRow(bvhNode, i);
    }
    return result;
}
function toBVH(data, fps = 24) {
    const bvhNode = subTree(data.joints);
    fillKeyFrames(data, bvhNode, fps);
    return bvhNode;
}


/***/ }),

/***/ "./src/hierarchy.ts":
/*!**************************!*\
  !*** ./src/hierarchy.ts ***!
  \**************************/
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.hierarchy = void 0;
exports.hierarchy = {
    "bvhName": "hip",
    "children": [
        {
            "bvhName": "mSpine1",
            "children": [
                {
                    "bvhName": "mSpine2",
                    "children": [
                        {
                            "bvhName": "abdomen",
                            "children": [
                                {
                                    "bvhName": "mSpine3",
                                    "children": [
                                        {
                                            "bvhName": "mSpine4",
                                            "children": [
                                                {
                                                    "bvhName": "chest",
                                                    "children": [
                                                        {
                                                            "bvhName": "neck",
                                                            "children": [
                                                                {
                                                                    "bvhName": "head",
                                                                    "children": [
                                                                        {
                                                                            "bvhName": "figureHair",
                                                                            "children": [
                                                                                {
                                                                                    "bvhName": "end"
                                                                                }
                                                                            ]
                                                                        },
                                                                        {
                                                                            "bvhName": "mEyeRight",
                                                                            "children": [
                                                                                {
                                                                                    "bvhName": "end"
                                                                                }
                                                                            ]
                                                                        },
                                                                        {
                                                                            "bvhName": "mEyeLeft",
                                                                            "children": [
                                                                                {
                                                                                    "bvhName": "end"
                                                                                }
                                                                            ]
                                                                        },
                                                                        {
                                                                            "bvhName": "mFaceRoot",
                                                                            "children": [
                                                                                {
                                                                                    "bvhName": "mFaceEyeAltRight",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceEyeAltLeft",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceForeheadLeft",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceForeheadRight",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceEyebrowOuterLeft",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceEyebrowCenterLeft",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceEyebrowInnerLeft",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceEyebrowOuterRight",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceEyebrowCenterRight",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceEyebrowInnerRight",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceEyeLidUpperLeft",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceEyeLidLowerLeft",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceEyeLidUpperRight",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceEyeLidLowerRight",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceEar1Left",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "mFaceEar2Left",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "end"
                                                                                                }
                                                                                            ]
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceEar1Right",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "mFaceEar2Right",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "end"
                                                                                                }
                                                                                            ]
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceNoseLeft",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceNoseCenter",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceNoseRight",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceCheekLowerLeft",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceCheekUpperLeft",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceCheekLowerRight",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceCheekUpperRight",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceJaw",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "mFaceChin",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "end"
                                                                                                }
                                                                                            ]
                                                                                        },
                                                                                        {
                                                                                            "bvhName": "mFaceTeethLower",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "mFaceLipLowerLeft",
                                                                                                    "children": [
                                                                                                        {
                                                                                                            "bvhName": "end"
                                                                                                        }
                                                                                                    ]
                                                                                                },
                                                                                                {
                                                                                                    "bvhName": "mFaceLipLowerRight",
                                                                                                    "children": [
                                                                                                        {
                                                                                                            "bvhName": "end"
                                                                                                        }
                                                                                                    ]
                                                                                                },
                                                                                                {
                                                                                                    "bvhName": "mFaceLipLowerCenter",
                                                                                                    "children": [
                                                                                                        {
                                                                                                            "bvhName": "end"
                                                                                                        }
                                                                                                    ]
                                                                                                },
                                                                                                {
                                                                                                    "bvhName": "mFaceTongueBase",
                                                                                                    "children": [
                                                                                                        {
                                                                                                            "bvhName": "mFaceTongueTip",
                                                                                                            "children": [
                                                                                                                {
                                                                                                                    "bvhName": "end"
                                                                                                                }
                                                                                                            ]
                                                                                                        }
                                                                                                    ]
                                                                                                }
                                                                                            ]
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceJawShaper",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceForeheadCenter",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceNoseBase",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceTeethUpper",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "mFaceLipUpperLeft",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "end"
                                                                                                }
                                                                                            ]
                                                                                        },
                                                                                        {
                                                                                            "bvhName": "mFaceLipUpperRight",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "end"
                                                                                                }
                                                                                            ]
                                                                                        },
                                                                                        {
                                                                                            "bvhName": "mFaceLipCornerLeft",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "end"
                                                                                                }
                                                                                            ]
                                                                                        },
                                                                                        {
                                                                                            "bvhName": "mFaceLipCornerRight",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "end"
                                                                                                }
                                                                                            ]
                                                                                        },
                                                                                        {
                                                                                            "bvhName": "mFaceLipUpperCenter",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "end"
                                                                                                }
                                                                                            ]
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceEyecornerInnerLeft",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceEyecornerInnerRight",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                },
                                                                                {
                                                                                    "bvhName": "mFaceNoseBridge",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "end"
                                                                                        }
                                                                                    ]
                                                                                }
                                                                            ]
                                                                        }
                                                                    ]
                                                                }
                                                            ]
                                                        },
                                                        {
                                                            "bvhName": "lCollar",
                                                            "children": [
                                                                {
                                                                    "bvhName": "lShldr",
                                                                    "children": [
                                                                        {
                                                                            "bvhName": "lForeArm",
                                                                            "children": [
                                                                                {
                                                                                    "bvhName": "lHand",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "mHandMiddle1Left",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "mHandMiddle2Left",
                                                                                                    "children": [
                                                                                                        {
                                                                                                            "bvhName": "mHandMiddle3Left",
                                                                                                            "children": [
                                                                                                                {
                                                                                                                    "bvhName": "end"
                                                                                                                }
                                                                                                            ]
                                                                                                        }
                                                                                                    ]
                                                                                                }
                                                                                            ]
                                                                                        },
                                                                                        {
                                                                                            "bvhName": "mHandIndex1Left",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "mHandIndex2Left",
                                                                                                    "children": [
                                                                                                        {
                                                                                                            "bvhName": "mHandIndex3Left",
                                                                                                            "children": [
                                                                                                                {
                                                                                                                    "bvhName": "end"
                                                                                                                }
                                                                                                            ]
                                                                                                        }
                                                                                                    ]
                                                                                                }
                                                                                            ]
                                                                                        },
                                                                                        {
                                                                                            "bvhName": "mHandRing1Left",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "mHandRing2Left",
                                                                                                    "children": [
                                                                                                        {
                                                                                                            "bvhName": "mHandRing3Left",
                                                                                                            "children": [
                                                                                                                {
                                                                                                                    "bvhName": "end"
                                                                                                                }
                                                                                                            ]
                                                                                                        }
                                                                                                    ]
                                                                                                }
                                                                                            ]
                                                                                        },
                                                                                        {
                                                                                            "bvhName": "mHandPinky1Left",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "mHandPinky2Left",
                                                                                                    "children": [
                                                                                                        {
                                                                                                            "bvhName": "mHandPinky3Left",
                                                                                                            "children": [
                                                                                                                {
                                                                                                                    "bvhName": "end"
                                                                                                                }
                                                                                                            ]
                                                                                                        }
                                                                                                    ]
                                                                                                }
                                                                                            ]
                                                                                        },
                                                                                        {
                                                                                            "bvhName": "mHandThumb1Left",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "mHandThumb2Left",
                                                                                                    "children": [
                                                                                                        {
                                                                                                            "bvhName": "mHandThumb3Left",
                                                                                                            "children": [
                                                                                                                {
                                                                                                                    "bvhName": "end"
                                                                                                                }
                                                                                                            ]
                                                                                                        }
                                                                                                    ]
                                                                                                }
                                                                                            ]
                                                                                        }
                                                                                    ]
                                                                                }
                                                                            ]
                                                                        }
                                                                    ]
                                                                }
                                                            ]
                                                        },
                                                        {
                                                            "bvhName": "rCollar",
                                                            "children": [
                                                                {
                                                                    "bvhName": "rShldr",
                                                                    "children": [
                                                                        {
                                                                            "bvhName": "rForeArm",
                                                                            "children": [
                                                                                {
                                                                                    "bvhName": "rHand",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "mHandMiddle1Right",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "mHandMiddle2Right",
                                                                                                    "children": [
                                                                                                        {
                                                                                                            "bvhName": "mHandMiddle3Right",
                                                                                                            "children": [
                                                                                                                {
                                                                                                                    "bvhName": "end"
                                                                                                                }
                                                                                                            ]
                                                                                                        }
                                                                                                    ]
                                                                                                }
                                                                                            ]
                                                                                        },
                                                                                        {
                                                                                            "bvhName": "mHandIndex1Right",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "mHandIndex2Right",
                                                                                                    "children": [
                                                                                                        {
                                                                                                            "bvhName": "mHandIndex3Right",
                                                                                                            "children": [
                                                                                                                {
                                                                                                                    "bvhName": "end"
                                                                                                                }
                                                                                                            ]
                                                                                                        }
                                                                                                    ]
                                                                                                }
                                                                                            ]
                                                                                        },
                                                                                        {
                                                                                            "bvhName": "mHandRing1Right",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "mHandRing2Right",
                                                                                                    "children": [
                                                                                                        {
                                                                                                            "bvhName": "mHandRing3Right",
                                                                                                            "children": [
                                                                                                                {
                                                                                                                    "bvhName": "end"
                                                                                                                }
                                                                                                            ]
                                                                                                        }
                                                                                                    ]
                                                                                                }
                                                                                            ]
                                                                                        },
                                                                                        {
                                                                                            "bvhName": "mHandPinky1Right",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "mHandPinky2Right",
                                                                                                    "children": [
                                                                                                        {
                                                                                                            "bvhName": "mHandPinky3Right",
                                                                                                            "children": [
                                                                                                                {
                                                                                                                    "bvhName": "end"
                                                                                                                }
                                                                                                            ]
                                                                                                        }
                                                                                                    ]
                                                                                                }
                                                                                            ]
                                                                                        },
                                                                                        {
                                                                                            "bvhName": "mHandThumb1Right",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "mHandThumb2Right",
                                                                                                    "children": [
                                                                                                        {
                                                                                                            "bvhName": "mHandThumb3Right",
                                                                                                            "children": [
                                                                                                                {
                                                                                                                    "bvhName": "end"
                                                                                                                }
                                                                                                            ]
                                                                                                        }
                                                                                                    ]
                                                                                                }
                                                                                            ]
                                                                                        }
                                                                                    ]
                                                                                }
                                                                            ]
                                                                        }
                                                                    ]
                                                                }
                                                            ]
                                                        },
                                                        {
                                                            "bvhName": "mWingsRoot",
                                                            "children": [
                                                                {
                                                                    "bvhName": "mWing1Left",
                                                                    "children": [
                                                                        {
                                                                            "bvhName": "mWing2Left",
                                                                            "children": [
                                                                                {
                                                                                    "bvhName": "mWing3Left",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "mWing4Left",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "end"
                                                                                                }
                                                                                            ]
                                                                                        },
                                                                                        {
                                                                                            "bvhName": "mWing4FanLeft",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "end"
                                                                                                }
                                                                                            ]
                                                                                        }
                                                                                    ]
                                                                                }
                                                                            ]
                                                                        }
                                                                    ]
                                                                },
                                                                {
                                                                    "bvhName": "mWing1Right",
                                                                    "children": [
                                                                        {
                                                                            "bvhName": "mWing2Right",
                                                                            "children": [
                                                                                {
                                                                                    "bvhName": "mWing3Right",
                                                                                    "children": [
                                                                                        {
                                                                                            "bvhName": "mWing4Right",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "end"
                                                                                                }
                                                                                            ]
                                                                                        },
                                                                                        {
                                                                                            "bvhName": "mWing4FanRight",
                                                                                            "children": [
                                                                                                {
                                                                                                    "bvhName": "end"
                                                                                                }
                                                                                            ]
                                                                                        }
                                                                                    ]
                                                                                }
                                                                            ]
                                                                        }
                                                                    ]
                                                                }
                                                            ]
                                                        }
                                                    ]
                                                }
                                            ]
                                        }
                                    ]
                                }
                            ]
                        }
                    ]
                }
            ]
        },
        {
            "bvhName": "rThigh",
            "children": [
                {
                    "bvhName": "rShin",
                    "children": [
                        {
                            "bvhName": "rFoot",
                            "children": [
                                {
                                    "bvhName": "mFootRight",
                                    "children": [
                                        {
                                            "bvhName": "mToeRight",
                                            "children": [
                                                {
                                                    "bvhName": "end"
                                                }
                                            ]
                                        }
                                    ]
                                }
                            ]
                        }
                    ]
                }
            ]
        },
        {
            "bvhName": "lThigh",
            "children": [
                {
                    "bvhName": "lShin",
                    "children": [
                        {
                            "bvhName": "lFoot",
                            "children": [
                                {
                                    "bvhName": "mFootLeft",
                                    "children": [
                                        {
                                            "bvhName": "mToeLeft",
                                            "children": [
                                                {
                                                    "bvhName": "end"
                                                }
                                            ]
                                        }
                                    ]
                                }
                            ]
                        }
                    ]
                }
            ]
        },
        {
            "bvhName": "mTail1",
            "children": [
                {
                    "bvhName": "mTail2",
                    "children": [
                        {
                            "bvhName": "mTail3",
                            "children": [
                                {
                                    "bvhName": "mTail4",
                                    "children": [
                                        {
                                            "bvhName": "mTail5",
                                            "children": [
                                                {
                                                    "bvhName": "mTail6",
                                                    "children": [
                                                        {
                                                            "bvhName": "end"
                                                        }
                                                    ]
                                                }
                                            ]
                                        }
                                    ]
                                }
                            ]
                        }
                    ]
                }
            ]
        },
        {
            "bvhName": "mGroin",
            "children": [
                {
                    "bvhName": "end"
                }
            ]
        },
        {
            "bvhName": "mHindLimbsRoot",
            "children": [
                {
                    "bvhName": "mHindLimb1Left",
                    "children": [
                        {
                            "bvhName": "mHindLimb2Left",
                            "children": [
                                {
                                    "bvhName": "mHindLimb3Left",
                                    "children": [
                                        {
                                            "bvhName": "mHindLimb4Left",
                                            "children": [
                                                {
                                                    "bvhName": "end"
                                                }
                                            ]
                                        }
                                    ]
                                }
                            ]
                        }
                    ]
                },
                {
                    "bvhName": "mHindLimb1Right",
                    "children": [
                        {
                            "bvhName": "mHindLimb2Right",
                            "children": [
                                {
                                    "bvhName": "mHindLimb3Right",
                                    "children": [
                                        {
                                            "bvhName": "mHindLimb4Right",
                                            "children": [
                                                {
                                                    "bvhName": "end"
                                                }
                                            ]
                                        }
                                    ]
                                }
                            ]
                        }
                    ]
                }
            ]
        }
    ]
};


/***/ }),

/***/ "./src/offsets.ts":
/*!************************!*\
  !*** ./src/offsets.ts ***!
  \************************/
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.maleOffsets = exports.femaleOffsets = void 0;
exports.femaleOffsets = {
    "hip": {
        "x": 0,
        "y": 0,
        "z": 0
    },
    "mHindLimbsRoot": {
        "x": 0,
        "y": 3.3071,
        "z": -7.874
    },
    "mHindLimb1Right": {
        "x": -5.0787,
        "y": -4.9213,
        "z": -8.0315
    },
    "mHindLimb2Right": {
        "x": 1.811,
        "y": -19.3307,
        "z": 0.0787
    },
    "mHindLimb3Right": {
        "x": 0.1181,
        "y": -18.4252,
        "z": -1.1811
    },
    "mHindLimb4Right": {
        "x": 0,
        "y": -2.4016,
        "z": 4.4094
    },
    "end_mHindLimb4Right": {
        "x": -0.315,
        "y": 0,
        "z": 4.1339
    },
    "mHindLimb1Left": {
        "x": 5.0787,
        "y": -4.9213,
        "z": -8.0315
    },
    "mHindLimb2Left": {
        "x": -1.811,
        "y": -19.3307,
        "z": 0.0787
    },
    "mHindLimb3Left": {
        "x": -0.1181,
        "y": -18.4252,
        "z": -1.1811
    },
    "mHindLimb4Left": {
        "x": 0,
        "y": -2.4016,
        "z": 4.4094
    },
    "end_mHindLimb4Left": {
        "x": 0.315,
        "y": 0,
        "z": 4.1339
    },
    "mGroin": {
        "x": 0,
        "y": -3.8189,
        "z": 2.5197
    },
    "end_mGroin": {
        "x": 0,
        "y": -2.5984,
        "z": 0.1575
    },
    "mTail1": {
        "x": 0,
        "y": 1.8504,
        "z": -4.5669
    },
    "mTail2": {
        "x": 0,
        "y": 0,
        "z": -7.7559
    },
    "mTail3": {
        "x": 0,
        "y": 0,
        "z": -6.6142
    },
    "mTail4": {
        "x": 0,
        "y": 0,
        "z": -5.5905
    },
    "mTail5": {
        "x": 0,
        "y": 0,
        "z": -4.4094
    },
    "mTail6": {
        "x": 0,
        "y": 0,
        "z": -3.7008
    },
    "end_mTail6": {
        "x": 0,
        "y": 0,
        "z": -3.5039
    },
    "lThigh": {
        "x": 4.9907,
        "y": -1.6141,
        "z": 1.329
    },
    "lShin": {
        "x": -1.794,
        "y": -19.3328,
        "z": -0.0349
    },
    "lFoot": {
        "x": 0.0543,
        "y": -18.4428,
        "z": -1.1373
    },
    "mFootLeft": {
        "x": 0,
        "y": -2.3866,
        "z": 4.4077
    },
    "mToeLeft": {
        "x": 0,
        "y": 0,
        "z": 4.2913
    },
    "end_mToeLeft": {
        "x": 0,
        "y": 0,
        "z": 0.7874
    },
    "rThigh": {
        "x": -5.0711,
        "y": -1.6176,
        "z": 1.3236
    },
    "rShin": {
        "x": 1.9148,
        "y": -19.3276,
        "z": -0.0307
    },
    "rFoot": {
        "x": 0,
        "y": -18.4446,
        "z": -1.1366
    },
    "mFootRight": {
        "x": 0,
        "y": -2.3873,
        "z": 4.4077
    },
    "mToeRight": {
        "x": 0,
        "y": 0,
        "z": 4.2913
    },
    "end_mToeRight": {
        "x": 0,
        "y": 0,
        "z": 0.7874
    },
    "mSpine1": {
        "x": 0,
        "y": 3.31,
        "z": 0
    },
    "mSpine2": {
        "x": 0,
        "y": -3.31,
        "z": 0
    },
    "abdomen": {
        "x": 0,
        "y": 3.31,
        "z": 0
    },
    "mSpine3": {
        "x": 0,
        "y": 8.066,
        "z": -0.605
    },
    "mSpine4": {
        "x": 0,
        "y": -8.066,
        "z": 0.605
    },
    "chest": {
        "x": 0,
        "y": 8.066,
        "z": -0.605
    },
    "mWingsRoot": {
        "x": 0,
        "y": 0,
        "z": -0.5457
    },
    "mWing1Right": {
        "x": -3.7402,
        "y": 7.126,
        "z": -3.8976
    },
    "mWing2Right": {
        "x": -6.6535,
        "y": 2.6378,
        "z": -6.6142
    },
    "mWing3Right": {
        "x": -7.2047,
        "y": 0,
        "z": -7.126
    },
    "mWing4FanRight": {
        "x": -6.811,
        "y": 0,
        "z": -6.7323
    },
    "end_mWing4FanRight": {
        "x": -2.4409,
        "y": -6.2598,
        "z": -2.6772
    },
    "mWing4Right": {
        "x": -6.811,
        "y": 0,
        "z": -6.7323
    },
    "end_mWing4Right": {
        "x": -5.1968,
        "y": 0,
        "z": -5.748
    },
    "mWing1Left": {
        "x": 3.7402,
        "y": 7.126,
        "z": -3.8976
    },
    "mWing2Left": {
        "x": 6.6535,
        "y": 2.6378,
        "z": -6.6142
    },
    "mWing3Left": {
        "x": 7.2047,
        "y": 0,
        "z": -7.126
    },
    "mWing4FanLeft": {
        "x": 6.811,
        "y": 0,
        "z": -6.7323
    },
    "end_mWing4FanLeft": {
        "x": 2.4409,
        "y": -6.2598,
        "z": -2.6772
    },
    "mWing4Left": {
        "x": 6.811,
        "y": 0,
        "z": -6.7323
    },
    "end_mWing4Left": {
        "x": 5.1968,
        "y": 0,
        "z": -5.748
    },
    "rCollar": {
        "x": -2.8346,
        "y": 6.5116,
        "z": -0.8157
    },
    "rShldr": {
        "x": -3.1267,
        "y": 0,
        "z": 0
    },
    "rForeArm": {
        "x": -10.9354,
        "y": 0,
        "z": 0
    },
    "rHand": {
        "x": -9.5236,
        "y": 0,
        "z": 0
    },
    "mHandThumb1Right": {
        "x": -1.0236,
        "y": 0.1575,
        "z": 1.2205
    },
    "mHandThumb2Right": {
        "x": -1.2598,
        "y": -0.0394,
        "z": 1.1024
    },
    "mHandThumb3Right": {
        "x": -1.2205,
        "y": -0.0394,
        "z": 0.9055
    },
    "end_mHandThumb3Right": {
        "x": -0.9843,
        "y": 0,
        "z": 0.5906
    },
    "mHandPinky1Right": {
        "x": -3.7401,
        "y": 0.1181,
        "z": -1.2205
    },
    "mHandPinky2Right": {
        "x": -0.9842,
        "y": -0.2362,
        "z": -0.9449
    },
    "mHandPinky3Right": {
        "x": -0.7087,
        "y": -0.1575,
        "z": -0.5906
    },
    "end_mHandPinky3Right": {
        "x": -0.6299,
        "y": -0.1575,
        "z": -0.5118
    },
    "mHandRing1Right": {
        "x": -3.8976,
        "y": 0.3543,
        "z": -0.3937
    },
    "mHandRing2Right": {
        "x": -1.4961,
        "y": -0.315,
        "z": -0.5118
    },
    "mHandRing3Right": {
        "x": -1.5748,
        "y": -0.3543,
        "z": -0.5118
    },
    "end_mHandRing3Right": {
        "x": -1.1024,
        "y": -0.2362,
        "z": -0.3937
    },
    "mHandIndex1Right": {
        "x": -3.8189,
        "y": 0.5905,
        "z": 1.4961
    },
    "mHandIndex2Right": {
        "x": -1.4173,
        "y": -0.2362,
        "z": 0.6693
    },
    "mHandIndex3Right": {
        "x": -1.2598,
        "y": -0.2362,
        "z": 0.5512
    },
    "end_mHandIndex3Right": {
        "x": -0.9843,
        "y": -0.1575,
        "z": 0.4331
    },
    "mHandMiddle1Right": {
        "x": -3.9764,
        "y": 0.5905,
        "z": 0.5118
    },
    "mHandMiddle2Right": {
        "x": -1.5748,
        "y": -0.2362,
        "z": -0.0394
    },
    "mHandMiddle3Right": {
        "x": -1.9291,
        "y": -0.315,
        "z": -0.0394
    },
    "end_mHandMiddle3Right": {
        "x": -1.2992,
        "y": -0.2362,
        "z": -0.0787
    },
    "lCollar": {
        "x": 2.822,
        "y": 6.5116,
        "z": -0.8157
    },
    "lShldr": {
        "x": 3.1102,
        "y": 0,
        "z": 0
    },
    "lForeArm": {
        "x": 10.9354,
        "y": 0,
        "z": 0
    },
    "lHand": {
        "x": 9.5165,
        "y": 0,
        "z": 0
    },
    "mHandThumb1Left": {
        "x": 1.0236,
        "y": 0.1575,
        "z": 1.2205
    },
    "mHandThumb2Left": {
        "x": 1.2598,
        "y": -0.0394,
        "z": 1.1024
    },
    "mHandThumb3Left": {
        "x": 1.2205,
        "y": -0.0394,
        "z": 0.9055
    },
    "end_mHandThumb3Left": {
        "x": 0.9843,
        "y": 0,
        "z": 0.5906
    },
    "mHandPinky1Left": {
        "x": 3.7402,
        "y": 0.1181,
        "z": -1.2205
    },
    "mHandPinky2Left": {
        "x": 0.9842,
        "y": -0.2362,
        "z": -0.9449
    },
    "mHandPinky3Left": {
        "x": 0.7087,
        "y": -0.1575,
        "z": -0.5906
    },
    "end_mHandPinky3Left": {
        "x": 0.6299,
        "y": -0.1575,
        "z": -0.5118
    },
    "mHandRing1Left": {
        "x": 3.8976,
        "y": 0.3543,
        "z": -0.3937
    },
    "mHandRing2Left": {
        "x": 1.4961,
        "y": -0.315,
        "z": -0.5118
    },
    "mHandRing3Left": {
        "x": 1.5748,
        "y": -0.3543,
        "z": -0.5118
    },
    "end_mHandRing3Left": {
        "x": 1.1024,
        "y": -0.2362,
        "z": -0.3937
    },
    "mHandIndex1Left": {
        "x": 3.8189,
        "y": 0.5905,
        "z": 1.4961
    },
    "mHandIndex2Left": {
        "x": 1.4173,
        "y": -0.2362,
        "z": 0.6693
    },
    "mHandIndex3Left": {
        "x": 1.2598,
        "y": -0.2362,
        "z": 0.5512
    },
    "end_mHandIndex3Left": {
        "x": 0.9843,
        "y": -0.1575,
        "z": 0.4331
    },
    "mHandMiddle1Left": {
        "x": 3.9764,
        "y": 0.5905,
        "z": 0.5118
    },
    "mHandMiddle2Left": {
        "x": 1.5748,
        "y": -0.2362,
        "z": -0.0394
    },
    "mHandMiddle3Left": {
        "x": 1.9291,
        "y": -0.315,
        "z": -0.0394
    },
    "end_mHandMiddle3Left": {
        "x": 1.2992,
        "y": -0.2362,
        "z": -0.0787
    },
    "neck": {
        "x": 0,
        "y": 9.8861,
        "z": -0.3705
    },
    "head": {
        "x": 0,
        "y": 2.9776,
        "z": 0
    },
    "mFaceRoot": {
        "x": 0,
        "y": 1.6934,
        "z": 0.9104
    },
    "mFaceNoseBridge": {
        "x": 0,
        "y": 0.7118,
        "z": 3.314
    },
    "end_mFaceNoseBridge": {
        "x": 0,
        "y": 0.315,
        "z": 0.5906
    },
    "mFaceEyecornerInnerRight": {
        "x": -0.7135,
        "y": 1.1122,
        "z": 2.7313
    },
    "end_mFaceEyecornerInnerRight": {
        "x": 0,
        "y": 0,
        "z": 0.6299
    },
    "mFaceEyecornerInnerLeft": {
        "x": 0.7135,
        "y": 1.1122,
        "z": 2.7313
    },
    "end_mFaceEyecornerInnerLeft": {
        "x": 0,
        "y": 0,
        "z": 0.6299
    },
    "mFaceTeethUpper": {
        "x": 0,
        "y": -1.1478,
        "z": 0.7283
    },
    "mFaceLipUpperCenter": {
        "x": 0,
        "y": -0.0964,
        "z": 1.6388
    },
    "end_mFaceLipUpperCenter": {
        "x": 0,
        "y": 0.0787,
        "z": 1.6929
    },
    "mFaceLipCornerRight": {
        "x": 0.6919,
        "y": -0.3642,
        "z": 1.0197
    },
    "end_mFaceLipCornerRight": {
        "x": -2.0079,
        "y": 0,
        "z": 1.7717
    },
    "mFaceLipCornerLeft": {
        "x": -0.6919,
        "y": -0.3642,
        "z": 1.0197
    },
    "end_mFaceLipCornerLeft": {
        "x": 2.0079,
        "y": 0,
        "z": 1.7717
    },
    "mFaceLipUpperRight": {
        "x": 0,
        "y": -0.1478,
        "z": 1.6388
    },
    "end_mFaceLipUpperRight": {
        "x": -0.5906,
        "y": 0,
        "z": 1.6142
    },
    "mFaceLipUpperLeft": {
        "x": 0,
        "y": -0.1478,
        "z": 1.6388
    },
    "end_mFaceLipUpperLeft": {
        "x": 0.5906,
        "y": 0,
        "z": 1.6142
    },
    "mFaceNoseBase": {
        "x": 0,
        "y": -0.7047,
        "z": 3.4232
    },
    "end_mFaceNoseBase": {
        "x": 0,
        "y": 0,
        "z": 0.5512
    },
    "mFaceForeheadCenter": {
        "x": 0,
        "y": 2.3027,
        "z": 2.5237
    },
    "end_mFaceForeheadCenter": {
        "x": 0,
        "y": 0,
        "z": 1.4173
    },
    "mFaceJawShaper": {
        "x": 0,
        "y": 0,
        "z": 0
    },
    "end_mFaceJawShaper": {
        "x": 0,
        "y": 0,
        "z": -0.6693
    },
    "mFaceJaw": {
        "x": 0,
        "y": -0.5339,
        "z": -0.0364
    },
    "mFaceTeethLower": {
        "x": 0,
        "y": -1.4869,
        "z": 0.7648
    },
    "mFaceTongueBase": {
        "x": 0,
        "y": 0.1821,
        "z": 1.4203
    },
    "mFaceTongueTip": {
        "x": 0,
        "y": 0.2549,
        "z": 0.8012
    },
    "end_mFaceTongueTip": {
        "x": 0,
        "y": 0,
        "z": 0.3937
    },
    "mFaceLipLowerCenter": {
        "x": 0,
        "y": 0,
        "z": 1.6388
    },
    "end_mFaceLipLowerCenter": {
        "x": 0,
        "y": 0.0787,
        "z": 1.5748
    },
    "mFaceLipLowerRight": {
        "x": 0,
        "y": 0,
        "z": 1.6388
    },
    "end_mFaceLipLowerRight": {
        "x": -0.6693,
        "y": 0.1969,
        "z": 1.3386
    },
    "mFaceLipLowerLeft": {
        "x": 0,
        "y": 0,
        "z": 1.6388
    },
    "end_mFaceLipLowerLeft": {
        "x": 0.6693,
        "y": 0.1969,
        "z": 1.3386
    },
    "mFaceChin": {
        "x": 0,
        "y": -2.0212,
        "z": 2.531
    },
    "end_mFaceChin": {
        "x": 0,
        "y": -0.7087,
        "z": 0.8268
    },
    "mFaceCheekUpperRight": {
        "x": -1.2663,
        "y": -0.2313,
        "z": 2.5492
    },
    "end_mFaceCheekUpperRight": {
        "x": -0.5906,
        "y": 0,
        "z": 0.8661
    },
    "mFaceCheekLowerRight": {
        "x": -1.2663,
        "y": -1.1941,
        "z": 1.8209
    },
    "end_mFaceCheekLowerRight": {
        "x": -1.1811,
        "y": 0,
        "z": 0.5118
    },
    "mFaceCheekUpperLeft": {
        "x": 1.2663,
        "y": -0.2313,
        "z": 2.5492
    },
    "end_mFaceCheekUpperLeft": {
        "x": 0.5906,
        "y": 0,
        "z": 0.8661
    },
    "mFaceCheekLowerLeft": {
        "x": 1.2663,
        "y": -1.1941,
        "z": 1.8209
    },
    "end_mFaceCheekLowerLeft": {
        "x": 1.1811,
        "y": 0,
        "z": 0.5118
    },
    "mFaceNoseRight": {
        "x": -0.5587,
        "y": -0.1957,
        "z": 3.1319
    },
    "end_mFaceNoseRight": {
        "x": -0.1575,
        "y": 0,
        "z": 0.5906
    },
    "mFaceNoseCenter": {
        "x": 0,
        "y": -0.0534,
        "z": 3.7146
    },
    "end_mFaceNoseCenter": {
        "x": 0,
        "y": 0,
        "z": 0.9843
    },
    "mFaceNoseLeft": {
        "x": 0.5587,
        "y": -0.1957,
        "z": 3.1319
    },
    "end_mFaceNoseLeft": {
        "x": 0.1575,
        "y": 0,
        "z": 0.5906
    },
    "mFaceEar1Right": {
        "x": -2.9795,
        "y": 0.0712,
        "z": 0
    },
    "mFaceEar2Right": {
        "x": -0.7087,
        "y": 0.9842,
        "z": -0.748
    },
    "end_mFaceEar2Right": {
        "x": 0,
        "y": 1.2992,
        "z": 0
    },
    "mFaceEar1Left": {
        "x": 2.9795,
        "y": 0.0712,
        "z": 0
    },
    "mFaceEar2Left": {
        "x": 0.7087,
        "y": 0.9842,
        "z": -0.748
    },
    "end_mFaceEar2Left": {
        "x": 0,
        "y": 1.2992,
        "z": 0
    },
    "mFaceEyeLidLowerRight": {
        "x": -1.4156,
        "y": 1.1834,
        "z": 2.6585
    },
    "end_mFaceEyeLidLowerRight": {
        "x": 0,
        "y": -0.2756,
        "z": 0.9449
    },
    "mFaceEyeLidUpperRight": {
        "x": -1.4156,
        "y": 1.1887,
        "z": 2.6585
    },
    "end_mFaceEyeLidUpperRight": {
        "x": 0,
        "y": 0.1969,
        "z": 1.063
    },
    "mFaceEyeLidLowerLeft": {
        "x": 1.4156,
        "y": 1.1834,
        "z": 2.6585
    },
    "end_mFaceEyeLidLowerLeft": {
        "x": 0,
        "y": -0.2756,
        "z": 0.9449
    },
    "mFaceEyeLidUpperLeft": {
        "x": 1.4156,
        "y": 1.1887,
        "z": 2.6585
    },
    "end_mFaceEyeLidUpperLeft": {
        "x": 0,
        "y": 0.1969,
        "z": 1.063
    },
    "mFaceEyebrowInnerRight": {
        "x": -0.906,
        "y": 1.8578,
        "z": 2.7313
    },
    "end_mFaceEyebrowInnerRight": {
        "x": 0,
        "y": 0,
        "z": 1.0236
    },
    "mFaceEyebrowCenterRight": {
        "x": -1.6711,
        "y": 2.1467,
        "z": 2.5492
    },
    "end_mFaceEyebrowCenterRight": {
        "x": 0,
        "y": 0,
        "z": 1.063
    },
    "mFaceEyebrowOuterRight": {
        "x": -1.9168,
        "y": 1.729,
        "z": 2.3307
    },
    "end_mFaceEyebrowOuterRight": {
        "x": -0.5118,
        "y": 0,
        "z": 0.9055
    },
    "mFaceEyebrowInnerLeft": {
        "x": 0.906,
        "y": 1.8578,
        "z": 2.7313
    },
    "end_mFaceEyebrowInnerLeft": {
        "x": 0,
        "y": 0,
        "z": 1.0236
    },
    "mFaceEyebrowCenterLeft": {
        "x": 1.6711,
        "y": 2.1467,
        "z": 2.5492
    },
    "end_mFaceEyebrowCenterLeft": {
        "x": 0,
        "y": 0,
        "z": 1.063
    },
    "mFaceEyebrowOuterLeft": {
        "x": 1.9168,
        "y": 1.729,
        "z": 2.3307
    },
    "end_mFaceEyebrowOuterLeft": {
        "x": 0.5118,
        "y": 0,
        "z": 0.9055
    },
    "mFaceForeheadRight": {
        "x": -1.3035,
        "y": 3.0608,
        "z": 2.3307
    },
    "end_mFaceForeheadRight": {
        "x": -0.1575,
        "y": 0.7087,
        "z": 0.9449
    },
    "mFaceForeheadLeft": {
        "x": 1.3035,
        "y": 3.0608,
        "z": 2.3307
    },
    "end_mFaceForeheadLeft": {
        "x": 0.1575,
        "y": 0.7087,
        "z": 0.9449
    },
    "mFaceEyeAltLeft": {
        "x": 1.4156,
        "y": 1.1836,
        "z": 2.6752
    },
    "end_mFaceEyeAltLeft": {
        "x": 0,
        "y": 0,
        "z": 0.9843
    },
    "mFaceEyeAltRight": {
        "x": -1.4156,
        "y": 1.1836,
        "z": 2.6754
    },
    "end_mFaceEyeAltRight": {
        "x": 0,
        "y": 0,
        "z": 0.9843
    },
    "mEyeLeft": {
        "x": 1.4203,
        "y": 2.877,
        "z": 3.5857
    },
    "end_mEyeLeft": {
        "x": 0,
        "y": 0,
        "z": 0.9843
    },
    "mEyeRight": {
        "x": -1.4203,
        "y": 2.877,
        "z": 3.5859
    },
    "end_mEyeRight": {
        "x": 0,
        "y": 0,
        "z": 0.9843
    },
    "figureHair": {
        "x": 0,
        "y": 2.6038,
        "z": 0
    },
    "end_figureHair": {
        "x": 0,
        "y": 1.2992,
        "z": 0
    }
};
exports.maleOffsets = {
    "hip": {
        "x": 0,
        "y": 0,
        "z": 0
    },
    "mHindLimbsRoot": {
        "x": 0,
        "y": 3.3071,
        "z": -7.874
    },
    "mHindLimb1Right": {
        "x": -5.0787,
        "y": -4.9213,
        "z": -8.0315
    },
    "mHindLimb2Right": {
        "x": 1.9016,
        "y": -19.3307,
        "z": 0.0827
    },
    "mHindLimb3Right": {
        "x": 0.124,
        "y": -20.2677,
        "z": -1.2402
    },
    "mHindLimb4Right": {
        "x": 0,
        "y": -2.4016,
        "z": 4.4094
    },
    "end_mHindLimb4Right": {
        "x": -0.315,
        "y": 0,
        "z": 4.1339
    },
    "mHindLimb1Left": {
        "x": 5.0787,
        "y": -4.9213,
        "z": -8.0315
    },
    "mHindLimb2Left": {
        "x": -1.9016,
        "y": -19.3307,
        "z": 0.0827
    },
    "mHindLimb3Left": {
        "x": -0.124,
        "y": -20.2677,
        "z": -1.2402
    },
    "mHindLimb4Left": {
        "x": 0,
        "y": -2.4016,
        "z": 4.4094
    },
    "end_mHindLimb4Left": {
        "x": 0.315,
        "y": 0,
        "z": 4.1339
    },
    "mGroin": {
        "x": 0,
        "y": -3.8189,
        "z": 2.5197
    },
    "end_mGroin": {
        "x": 0,
        "y": -2.5984,
        "z": 0.1575
    },
    "mTail1": {
        "x": 0,
        "y": 1.8504,
        "z": -4.5669
    },
    "mTail2": {
        "x": 0,
        "y": 0,
        "z": -7.7559
    },
    "mTail3": {
        "x": 0,
        "y": 0,
        "z": -6.6142
    },
    "mTail4": {
        "x": 0,
        "y": 0,
        "z": -5.5905
    },
    "mTail5": {
        "x": 0,
        "y": 0,
        "z": -4.4094
    },
    "mTail6": {
        "x": 0,
        "y": 0,
        "z": -3.7008
    },
    "end_mTail6": {
        "x": 0,
        "y": 0,
        "z": -3.5039
    },
    "lThigh": {
        "x": 4.9907,
        "y": -1.6141,
        "z": 1.329
    },
    "lShin": {
        "x": -1.8837,
        "y": -19.3328,
        "z": -0.0367
    },
    "lFoot": {
        "x": 0.057,
        "y": -20.2871,
        "z": -1.1941
    },
    "mFootLeft": {
        "x": 0,
        "y": -2.3866,
        "z": 4.4077
    },
    "mToeLeft": {
        "x": 0,
        "y": 0,
        "z": 4.2913
    },
    "end_mToeLeft": {
        "x": 0,
        "y": 0,
        "z": 0.7874
    },
    "rThigh": {
        "x": -5.0711,
        "y": -1.6176,
        "z": 1.3236
    },
    "rShin": {
        "x": 2.0105,
        "y": -19.3276,
        "z": -0.0322
    },
    "rFoot": {
        "x": 0,
        "y": -20.2891,
        "z": -1.1934
    },
    "mFootRight": {
        "x": 0,
        "y": -2.3873,
        "z": 4.4077
    },
    "mToeRight": {
        "x": 0,
        "y": 0,
        "z": 4.2913
    },
    "end_mToeRight": {
        "x": 0,
        "y": 0,
        "z": 0.7874
    },
    "mSpine1": {
        "x": 0,
        "y": 3.31,
        "z": 0
    },
    "mSpine2": {
        "x": 0,
        "y": -3.31,
        "z": 0
    },
    "abdomen": {
        "x": 0,
        "y": 3.31,
        "z": 0
    },
    "mSpine3": {
        "x": 0,
        "y": 8.4693,
        "z": -0.605
    },
    "mSpine4": {
        "x": 0,
        "y": -8.4693,
        "z": 0.605
    },
    "chest": {
        "x": 0,
        "y": 8.4693,
        "z": -0.605
    },
    "mWingsRoot": {
        "x": 0,
        "y": 0,
        "z": -0.5732
    },
    "mWing1Right": {
        "x": -3.7402,
        "y": 7.126,
        "z": -5.8661
    },
    "mWing2Right": {
        "x": -6.6535,
        "y": 2.6378,
        "z": -6.6142
    },
    "mWing3Right": {
        "x": -7.2047,
        "y": 0,
        "z": -7.126
    },
    "mWing4FanRight": {
        "x": -6.811,
        "y": 0,
        "z": -6.7323
    },
    "end_mWing4FanRight": {
        "x": -2.4409,
        "y": -6.2598,
        "z": -2.6772
    },
    "mWing4Right": {
        "x": -6.811,
        "y": 0,
        "z": -6.7323
    },
    "end_mWing4Right": {
        "x": -5.1968,
        "y": 0,
        "z": -5.748
    },
    "mWing1Left": {
        "x": 3.7402,
        "y": 7.126,
        "z": -5.8661
    },
    "mWing2Left": {
        "x": 6.6535,
        "y": 2.6378,
        "z": -6.6142
    },
    "mWing3Left": {
        "x": 7.2047,
        "y": 0,
        "z": -7.126
    },
    "mWing4FanLeft": {
        "x": 6.811,
        "y": 0,
        "z": -6.7323
    },
    "end_mWing4FanLeft": {
        "x": 2.4409,
        "y": -6.2598,
        "z": -2.6772
    },
    "mWing4Left": {
        "x": 6.811,
        "y": 0,
        "z": -6.7323
    },
    "end_mWing4Left": {
        "x": 5.1968,
        "y": 0,
        "z": -5.748
    },
    "rCollar": {
        "x": -2.9823,
        "y": 6.8372,
        "z": -0.8569
    },
    "rShldr": {
        "x": -4.3774,
        "y": 0,
        "z": 0
    },
    "rForeArm": {
        "x": -14.3527,
        "y": 0,
        "z": 0
    },
    "rHand": {
        "x": -10.3307,
        "y": 0,
        "z": 0
    },
    "mHandThumb1Right": {
        "x": -1.0236,
        "y": 0.1575,
        "z": 1.2205
    },
    "mHandThumb2Right": {
        "x": -1.2598,
        "y": -0.0394,
        "z": 1.1024
    },
    "mHandThumb3Right": {
        "x": -1.2205,
        "y": -0.0394,
        "z": 0.9055
    },
    "end_mHandThumb3Right": {
        "x": -0.9843,
        "y": 0,
        "z": 0.5906
    },
    "mHandPinky1Right": {
        "x": -3.7402,
        "y": 0.1181,
        "z": -1.2205
    },
    "mHandPinky2Right": {
        "x": -0.9842,
        "y": -0.2362,
        "z": -0.9449
    },
    "mHandPinky3Right": {
        "x": -0.7087,
        "y": -0.1575,
        "z": -0.5906
    },
    "end_mHandPinky3Right": {
        "x": -0.6299,
        "y": -0.1575,
        "z": -0.5118
    },
    "mHandRing1Right": {
        "x": -3.8976,
        "y": 0.3543,
        "z": -0.3937
    },
    "mHandRing2Right": {
        "x": -1.4961,
        "y": -0.315,
        "z": -0.5118
    },
    "mHandRing3Right": {
        "x": -1.5748,
        "y": -0.3543,
        "z": -0.5118
    },
    "end_mHandRing3Right": {
        "x": -1.1024,
        "y": -0.2362,
        "z": -0.3937
    },
    "mHandIndex1Right": {
        "x": -3.8189,
        "y": 0.5905,
        "z": 1.4961
    },
    "mHandIndex2Right": {
        "x": -1.4173,
        "y": -0.2362,
        "z": 0.6693
    },
    "mHandIndex3Right": {
        "x": -1.2598,
        "y": -0.2362,
        "z": 0.5512
    },
    "end_mHandIndex3Right": {
        "x": -0.9843,
        "y": -0.1575,
        "z": 0.4331
    },
    "mHandMiddle1Right": {
        "x": -3.9764,
        "y": 0.5905,
        "z": 0.5118
    },
    "mHandMiddle2Right": {
        "x": -1.5748,
        "y": -0.2362,
        "z": -0.0394
    },
    "mHandMiddle3Right": {
        "x": -1.9291,
        "y": -0.315,
        "z": -0.0394
    },
    "end_mHandMiddle3Right": {
        "x": -1.2992,
        "y": -0.2362,
        "z": -0.0787
    },
    "lCollar": {
        "x": 2.969,
        "y": 6.8372,
        "z": -0.8569
    },
    "lShldr": {
        "x": 4.3543,
        "y": 0,
        "z": 0
    },
    "lForeArm": {
        "x": 14.3527,
        "y": 0,
        "z": 0
    },
    "lHand": {
        "x": 10.3229,
        "y": 0,
        "z": 0
    },
    "mHandThumb1Left": {
        "x": 1.0236,
        "y": 0.1575,
        "z": 1.2205
    },
    "mHandThumb2Left": {
        "x": 1.2598,
        "y": -0.0394,
        "z": 1.1024
    },
    "mHandThumb3Left": {
        "x": 1.2205,
        "y": -0.0394,
        "z": 0.9055
    },
    "end_mHandThumb3Left": {
        "x": 0.9843,
        "y": 0,
        "z": 0.5906
    },
    "mHandPinky1Left": {
        "x": 3.7402,
        "y": 0.1181,
        "z": -1.2205
    },
    "mHandPinky2Left": {
        "x": 0.9842,
        "y": -0.2362,
        "z": -0.9449
    },
    "mHandPinky3Left": {
        "x": 0.7087,
        "y": -0.1575,
        "z": -0.5906
    },
    "end_mHandPinky3Left": {
        "x": 0.6299,
        "y": -0.1575,
        "z": -0.5118
    },
    "mHandRing1Left": {
        "x": 3.8976,
        "y": 0.3543,
        "z": -0.3937
    },
    "mHandRing2Left": {
        "x": 1.4961,
        "y": -0.315,
        "z": -0.5118
    },
    "mHandRing3Left": {
        "x": 1.5748,
        "y": -0.3543,
        "z": -0.5118
    },
    "end_mHandRing3Left": {
        "x": 1.1024,
        "y": -0.2362,
        "z": -0.3937
    },
    "mHandIndex1Left": {
        "x": 3.8189,
        "y": 0.5905,
        "z": 1.4961
    },
    "mHandIndex2Left": {
        "x": 1.4173,
        "y": -0.2362,
        "z": 0.6693
    },
    "mHandIndex3Left": {
        "x": 1.2598,
        "y": -0.2362,
        "z": 0.5512
    },
    "end_mHandIndex3Left": {
        "x": 0.9843,
        "y": -0.1575,
        "z": 0.4331
    },
    "mHandMiddle1Left": {
        "x": 3.9764,
        "y": 0.5905,
        "z": 0.5118
    },
    "mHandMiddle2Left": {
        "x": 1.5748,
        "y": -0.2362,
        "z": -0.0394
    },
    "mHandMiddle3Left": {
        "x": 1.9291,
        "y": -0.315,
        "z": -0.0394
    },
    "end_mHandMiddle3Left": {
        "x": 1.2992,
        "y": -0.2362,
        "z": -0.0787
    },
    "neck": {
        "x": 0,
        "y": 10.3804,
        "z": -0.3893
    },
    "head": {
        "x": 0,
        "y": 3.5731,
        "z": 0
    },
    "mFaceRoot": {
        "x": 0,
        "y": 1.6934,
        "z": 0.9104
    },
    "mFaceNoseBridge": {
        "x": 0,
        "y": 0.7118,
        "z": 3.314
    },
    "end_mFaceNoseBridge": {
        "x": 0,
        "y": 0.315,
        "z": 0.5906
    },
    "mFaceEyecornerInnerRight": {
        "x": -0.7135,
        "y": 1.1122,
        "z": 2.7313
    },
    "end_mFaceEyecornerInnerRight": {
        "x": 0,
        "y": 0,
        "z": 0.6299
    },
    "mFaceEyecornerInnerLeft": {
        "x": 0.7135,
        "y": 1.1122,
        "z": 2.7313
    },
    "end_mFaceEyecornerInnerLeft": {
        "x": 0,
        "y": 0,
        "z": 0.6299
    },
    "mFaceTeethUpper": {
        "x": 0,
        "y": -1.1478,
        "z": 0.7283
    },
    "mFaceLipUpperCenter": {
        "x": 0,
        "y": -0.0964,
        "z": 1.6388
    },
    "end_mFaceLipUpperCenter": {
        "x": 0,
        "y": 0.0787,
        "z": 1.6929
    },
    "mFaceLipCornerRight": {
        "x": 0.6919,
        "y": -0.3642,
        "z": 1.0197
    },
    "end_mFaceLipCornerRight": {
        "x": -2.0079,
        "y": 0,
        "z": 1.7717
    },
    "mFaceLipCornerLeft": {
        "x": -0.6919,
        "y": -0.3642,
        "z": 1.0197
    },
    "end_mFaceLipCornerLeft": {
        "x": 2.0079,
        "y": 0,
        "z": 1.7717
    },
    "mFaceLipUpperRight": {
        "x": 0,
        "y": -0.1478,
        "z": 1.6388
    },
    "end_mFaceLipUpperRight": {
        "x": -0.5906,
        "y": 0,
        "z": 1.6142
    },
    "mFaceLipUpperLeft": {
        "x": 0,
        "y": -0.1478,
        "z": 1.6388
    },
    "end_mFaceLipUpperLeft": {
        "x": 0.5906,
        "y": 0,
        "z": 1.6142
    },
    "mFaceNoseBase": {
        "x": 0,
        "y": -0.7047,
        "z": 3.4232
    },
    "end_mFaceNoseBase": {
        "x": 0,
        "y": 0,
        "z": 0.5512
    },
    "mFaceForeheadCenter": {
        "x": 0,
        "y": 2.3027,
        "z": 2.5237
    },
    "end_mFaceForeheadCenter": {
        "x": 0,
        "y": 0,
        "z": 1.4173
    },
    "mFaceJawShaper": {
        "x": 0,
        "y": 0,
        "z": 0
    },
    "end_mFaceJawShaper": {
        "x": 0,
        "y": 0,
        "z": -0.6693
    },
    "mFaceJaw": {
        "x": 0,
        "y": -0.5339,
        "z": -0.0364
    },
    "mFaceTeethLower": {
        "x": 0,
        "y": -1.4869,
        "z": 0.7648
    },
    "mFaceTongueBase": {
        "x": 0,
        "y": 0.1821,
        "z": 1.4203
    },
    "mFaceTongueTip": {
        "x": 0,
        "y": 0.2549,
        "z": 0.8012
    },
    "end_mFaceTongueTip": {
        "x": 0,
        "y": 0,
        "z": 0.3937
    },
    "mFaceLipLowerCenter": {
        "x": 0,
        "y": 0,
        "z": 1.6388
    },
    "end_mFaceLipLowerCenter": {
        "x": 0,
        "y": 0.0787,
        "z": 1.5748
    },
    "mFaceLipLowerRight": {
        "x": 0,
        "y": 0,
        "z": 1.6388
    },
    "end_mFaceLipLowerRight": {
        "x": -0.6693,
        "y": 0.1969,
        "z": 1.3386
    },
    "mFaceLipLowerLeft": {
        "x": 0,
        "y": 0,
        "z": 1.6388
    },
    "end_mFaceLipLowerLeft": {
        "x": 0.6693,
        "y": 0.1969,
        "z": 1.3386
    },
    "mFaceChin": {
        "x": 0,
        "y": -2.0212,
        "z": 2.531
    },
    "end_mFaceChin": {
        "x": 0,
        "y": -0.7087,
        "z": 0.8268
    },
    "mFaceCheekUpperRight": {
        "x": -1.2663,
        "y": -0.2313,
        "z": 2.5492
    },
    "end_mFaceCheekUpperRight": {
        "x": -0.5906,
        "y": 0,
        "z": 0.8661
    },
    "mFaceCheekLowerRight": {
        "x": -1.2663,
        "y": -1.1941,
        "z": 1.8209
    },
    "end_mFaceCheekLowerRight": {
        "x": -1.1811,
        "y": 0,
        "z": 0.5118
    },
    "mFaceCheekUpperLeft": {
        "x": 1.2663,
        "y": -0.2313,
        "z": 2.5492
    },
    "end_mFaceCheekUpperLeft": {
        "x": 0.5906,
        "y": 0,
        "z": 0.8661
    },
    "mFaceCheekLowerLeft": {
        "x": 1.2663,
        "y": -1.1941,
        "z": 1.8209
    },
    "end_mFaceCheekLowerLeft": {
        "x": 1.1811,
        "y": 0,
        "z": 0.5118
    },
    "mFaceNoseRight": {
        "x": -0.5587,
        "y": -0.1957,
        "z": 3.1319
    },
    "end_mFaceNoseRight": {
        "x": -0.1575,
        "y": 0,
        "z": 0.5906
    },
    "mFaceNoseCenter": {
        "x": 0,
        "y": -0.0534,
        "z": 3.7146
    },
    "end_mFaceNoseCenter": {
        "x": 0,
        "y": 0,
        "z": 0.9843
    },
    "mFaceNoseLeft": {
        "x": 0.5587,
        "y": -0.1957,
        "z": 3.1319
    },
    "end_mFaceNoseLeft": {
        "x": 0.1575,
        "y": 0,
        "z": 0.5906
    },
    "mFaceEar1Right": {
        "x": -2.9795,
        "y": 0.0712,
        "z": 0
    },
    "mFaceEar2Right": {
        "x": -0.7087,
        "y": 0.9842,
        "z": -0.748
    },
    "end_mFaceEar2Right": {
        "x": 0,
        "y": 1.2992,
        "z": 0
    },
    "mFaceEar1Left": {
        "x": 2.9795,
        "y": 0.0712,
        "z": 0
    },
    "mFaceEar2Left": {
        "x": 0.7087,
        "y": 0.9842,
        "z": -0.748
    },
    "end_mFaceEar2Left": {
        "x": 0,
        "y": 1.2992,
        "z": 0
    },
    "mFaceEyeLidLowerRight": {
        "x": -1.4156,
        "y": 1.1834,
        "z": 2.6585
    },
    "end_mFaceEyeLidLowerRight": {
        "x": 0,
        "y": -0.2756,
        "z": 0.9449
    },
    "mFaceEyeLidUpperRight": {
        "x": -1.4156,
        "y": 1.1887,
        "z": 2.6585
    },
    "end_mFaceEyeLidUpperRight": {
        "x": 0,
        "y": 0.1969,
        "z": 1.063
    },
    "mFaceEyeLidLowerLeft": {
        "x": 1.4156,
        "y": 1.1834,
        "z": 2.6585
    },
    "end_mFaceEyeLidLowerLeft": {
        "x": 0,
        "y": -0.2756,
        "z": 0.9449
    },
    "mFaceEyeLidUpperLeft": {
        "x": 1.4156,
        "y": 1.1887,
        "z": 2.6585
    },
    "end_mFaceEyeLidUpperLeft": {
        "x": 0,
        "y": 0.1969,
        "z": 1.063
    },
    "mFaceEyebrowInnerRight": {
        "x": -0.906,
        "y": 1.8578,
        "z": 2.7313
    },
    "end_mFaceEyebrowInnerRight": {
        "x": 0,
        "y": 0,
        "z": 1.0236
    },
    "mFaceEyebrowCenterRight": {
        "x": -1.6711,
        "y": 2.1467,
        "z": 2.5492
    },
    "end_mFaceEyebrowCenterRight": {
        "x": 0,
        "y": 0,
        "z": 1.063
    },
    "mFaceEyebrowOuterRight": {
        "x": -1.9168,
        "y": 1.729,
        "z": 2.3307
    },
    "end_mFaceEyebrowOuterRight": {
        "x": -0.5118,
        "y": 0,
        "z": 0.9055
    },
    "mFaceEyebrowInnerLeft": {
        "x": 0.906,
        "y": 1.8578,
        "z": 2.7313
    },
    "end_mFaceEyebrowInnerLeft": {
        "x": 0,
        "y": 0,
        "z": 1.0236
    },
    "mFaceEyebrowCenterLeft": {
        "x": 1.6711,
        "y": 2.1467,
        "z": 2.5492
    },
    "end_mFaceEyebrowCenterLeft": {
        "x": 0,
        "y": 0,
        "z": 1.063
    },
    "mFaceEyebrowOuterLeft": {
        "x": 1.9168,
        "y": 1.729,
        "z": 2.3307
    },
    "end_mFaceEyebrowOuterLeft": {
        "x": 0.5118,
        "y": 0,
        "z": 0.9055
    },
    "mFaceForeheadRight": {
        "x": -1.3035,
        "y": 3.0608,
        "z": 2.3307
    },
    "end_mFaceForeheadRight": {
        "x": -0.1575,
        "y": 0.7087,
        "z": 0.9449
    },
    "mFaceForeheadLeft": {
        "x": 1.3035,
        "y": 3.0608,
        "z": 2.3307
    },
    "end_mFaceForeheadLeft": {
        "x": 0.1575,
        "y": 0.7087,
        "z": 0.9449
    },
    "mFaceEyeAltLeft": {
        "x": 1.4156,
        "y": 1.1836,
        "z": 2.6752
    },
    "end_mFaceEyeAltLeft": {
        "x": 0,
        "y": 0,
        "z": 0.9843
    },
    "mFaceEyeAltRight": {
        "x": -1.4156,
        "y": 1.1836,
        "z": 2.6754
    },
    "end_mFaceEyeAltRight": {
        "x": 0,
        "y": 0,
        "z": 0.9843
    },
    "mEyeLeft": {
        "x": 1.4203,
        "y": 2.877,
        "z": 3.5857
    },
    "end_mEyeLeft": {
        "x": 0,
        "y": 0,
        "z": 0.9843
    },
    "mEyeRight": {
        "x": -1.4203,
        "y": 2.877,
        "z": 3.5859
    },
    "end_mEyeRight": {
        "x": 0,
        "y": 0,
        "z": 0.9843
    },
    "figureHair": {
        "x": 0,
        "y": 2.6038,
        "z": 0
    },
    "end_figureHair": {
        "x": 0,
        "y": 1.2992,
        "z": 0
    }
};


/***/ }),

/***/ "./src/parse.ts":
/*!**********************!*\
  !*** ./src/parse.ts ***!
  \**********************/
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.parseAnim = parseAnim;
exports.distributeSingleFrame = distributeSingleFrame;
exports.parseBVH = parseBVH;
const utils_1 = __webpack_require__(/*! ./utils */ "./src/utils.ts");
function parseAnim(arrayBuffer) {
    const view = new DataView(arrayBuffer);
    let offset = 0;
    function readU16() {
        const value = view.getUint16(offset, true);
        offset += 2;
        return value;
    }
    function readU32() {
        const value = view.getUint16(offset, true);
        offset += 4;
        return value;
    }
    function readS32() {
        const value = view.getInt32(offset, true);
        offset += 4;
        return value;
    }
    function readF32() {
        const value = view.getFloat32(offset, true);
        offset += 4;
        return value;
    }
    function readString() {
        let str = "";
        while (offset < view.byteLength) {
            const char = view.getUint8(offset++);
            if (char === 0)
                break; // NULL-terminated
            str += String.fromCharCode(char);
        }
        return str;
    }
    function readFixedString(size) {
        let str = "";
        for (let i = 0; i < size; i++) {
            const char = view.getUint8(offset++);
            if (char !== 0)
                str += String.fromCharCode(char);
        }
        return str;
    }
    const version = readU16();
    const sub_version = readU16();
    const base_priority = readS32();
    const duration = readF32();
    const emote_name = readString();
    const loop_in_point = readF32();
    const loop_out_point = readF32();
    const loop = readS32();
    const ease_in_duration = readF32();
    const ease_out_duration = readF32();
    const hand_pose = readU32();
    const num_joints = readU32();
    const joints = [];
    for (let i = 0; i < num_joints; i++) {
        const joint_name = readString();
        const joint_priority = readS32();
        const num_rot_keys = readS32();
        const rotation_keys = [];
        for (let j = 0; j < num_rot_keys; j++) {
            const time = readU16();
            const rot_x = readU16();
            const rot_y = readU16();
            const rot_z = readU16();
            rotation_keys.push({ time: intToFloat(time, 0, duration), x: intToFloat(rot_x, -1, 1), y: intToFloat(rot_y, -1, 1), z: intToFloat(rot_z, -1, 1) });
        }
        const num_pos_keys = readS32();
        const position_keys = [];
        for (let j = 0; j < num_pos_keys; j++) {
            const time = readU16();
            const pos_x = readU16();
            const pos_y = readU16();
            const pos_z = readU16();
            position_keys.push({ time: intToFloat(time, 0, duration), x: intToFloat(pos_x, -5, 5), y: intToFloat(pos_y, -5, 5), z: intToFloat(pos_z, -5, 5) });
        }
        joints.push({ joint_name, joint_priority, rotation_keys, position_keys });
    }
    const num_constraints = readS32();
    const constraints = [];
    for (let i = 0; i < num_constraints; i++) {
        const chain_length = view.getUint8(offset++);
        const constraint_type = view.getUint8(offset++);
        const source_volume = readFixedString(16);
        const source_offset = [readF32(), readF32(), readF32()];
        const target_volume = readFixedString(16);
        const target_offset = [readF32(), readF32(), readF32()];
        const target_dir = [readF32(), readF32(), readF32()];
        const ease_in_start = readF32();
        const ease_in_stop = readF32();
        const ease_out_start = readF32();
        const ease_out_stop = readF32();
        constraints.push({
            chain_length, constraint_type, source_volume, source_offset,
            target_volume, target_offset, target_dir, ease_in_start, ease_in_stop,
            ease_out_start, ease_out_stop
        });
    }
    joints.forEach((item) => item.rotation_keys.forEach((rot) => {
        if (!item.euler_keys) {
            item.euler_keys = [];
        }
        item.euler_keys.push((0, utils_1.toEulers)(rot));
    }));
    return { version, sub_version, duration, emote_name, loop, joints, constraints };
}
function intToFloat(val, min, max) {
    const one = (max - min) / 65535.0;
    const result = min + val * one;
    if (Math.abs(result) < one) {
        return 0;
    }
    return result;
}
function enumerate(content, key, alter) {
    let result = content;
    let count = 0;
    while (result.includes("\"" + key + "\"")) {
        result = result.replace("\"" + key + "\"", "\"" + alter + count + "\"");
        count += 1;
    }
    return result;
}
function compose(node, name) {
    const children = [];
    const result = {};
    let cnts = [];
    let jnts = [];
    if (Object.keys(node).includes("End Site")) {
        node.jnt1 = "end";
    }
    Object.keys(node).forEach(item => {
        if (item.includes("cnt")) {
            cnts.push(item);
            return;
        }
        if (item.includes("jnt")) {
            jnts.push(item);
        }
    });
    cnts = cnts.sort((val1, val2) => { return parseInt(val1.replace("cnt", "")) - parseInt(val2.replace("cnt", "")); });
    jnts = jnts.sort((val1, val2) => { return parseInt(val1.replace("jnt", "")) - parseInt(val2.replace("jnt", "")); });
    jnts.forEach((item, i) => children.push(compose(node[cnts[i]], node[item].trim())));
    if (node.OFFSET) {
        const offset = node.OFFSET.trim().split(" ").filter((item) => !isNaN(parseFloat(item))).map(parseFloat);
        result.offset = { x: offset[0], y: offset[1], z: offset[2] };
    }
    if (node.CHANNELS) {
        const channels = node.CHANNELS.trim().split(" ").map((item) => item.trim()).filter((item) => !!item).filter((item) => isNaN(parseFloat(item)));
        result.channels = channels;
    }
    result.bvhName = name;
    if (children.length) {
        result.children = children;
    }
    return result;
}
function cleanup(data) {
    data.jnt0 = data.ROOT;
    return compose(data, "root");
}
function parseFrames(rows) {
    const splitedRows = rows.map(item => item.split(" ").map(item => item.trim()).filter(item => !!item));
    return splitedRows.map(item => item.map(parseFloat));
}
function parseFramesPart(framesPart) {
    const framesRows = framesPart.split("\n");
    let timeIndex = -1;
    for (let i = 0; i < framesRows.length; i++) {
        if (framesRows[i].toLowerCase().includes("time")) {
            timeIndex = i;
            break;
        }
    }
    if (timeIndex < 0) {
        return { framesLength: 0, frameDuration: 0, frames: [] };
    }
    const framesLength = parseInt(framesRows[timeIndex - 1].split(" ").map((item) => item.trim()).filter((item) => !!item).filter((item) => !isNaN(item))[0]);
    const frameDuration = parseFloat(framesRows[timeIndex].split(" ").map((item) => item.trim()).filter((item) => !!item).filter((item) => !isNaN(item))[0]);
    while (!framesRows[0].toLowerCase().includes("time")) {
        framesRows.shift();
    }
    framesRows.shift();
    const frames = parseFrames(framesRows);
    return { framesLength, frameDuration, frames };
}
function distributeSingleFrame(hierarchy, frame) {
    var _a, _b;
    (_a = hierarchy.children) === null || _a === void 0 ? void 0 : _a.toReversed().forEach((child) => distributeSingleFrame(child, frame));
    const position = { x: 0, y: 0, z: 0 };
    const rotation = { x: 0, y: 0, z: 0 };
    (_b = hierarchy.channels) === null || _b === void 0 ? void 0 : _b.toReversed().forEach((item) => {
        const value = frame.pop();
        (0, utils_1.distributeValue)(position, rotation, item, value);
    });
    if (!hierarchy.bvhFrames) {
        hierarchy.bvhFrames = [];
    }
    hierarchy.bvhFrames.push({ position, rotation });
}
function parseBVH(text) {
    const parts = text.replaceAll("\r", "").replaceAll("\t", " ").split("MOTION");
    let result = parts[0].split("HIERARCHY")[1];
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789 ".split("").forEach(item => {
        result = result.replaceAll(item + "\n", item + "\",\n");
    });
    result = result.replaceAll("JOINT", "\"JOINT\":\"");
    result = result.replaceAll("OFFSET", "\"OFFSET\":\"");
    result = result.replaceAll("CHANNELS", "\"CHANNELS\":\"");
    result = result.replaceAll("ROOT", "\"ROOT\":\"");
    result = result.replaceAll("End Site", "\"End Site\":\"");
    result = result.replaceAll("{", "\"content\": {");
    result = result.split("}").map(item => {
        if (item.trim().endsWith(",")) {
            return item.trim() + "\"dummy\": {}";
        }
        return item;
    }).join("}");
    result = result.split("}").map(item => {
        if (item.trim().startsWith("\"JOINT\"")) {
            return "," + item.trim();
        }
        return item;
    }).join("}");
    let count = 0;
    result = enumerate(result, "JOINT", "jnt");
    result = enumerate(result, "content", "cnt");
    const hierarchy = cleanup(JSON.parse("{" + result + "}")).children[0];
    const animation = parseFramesPart(parts[1]);
    hierarchy.bvhTimes = [];
    for (let i = 0; i < animation.framesLength; i++) {
        hierarchy.bvhTimes.push(animation.frameDuration * i);
    }
    animation.frames.forEach(item => distributeSingleFrame(hierarchy, item));
    return hierarchy;
}


/***/ }),

/***/ "./src/utils.ts":
/*!**********************!*\
  !*** ./src/utils.ts ***!
  \**********************/
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.RAD_TO_DEG = void 0;
exports.toQuaternion = toQuaternion;
exports.toEulers = toEulers;
exports.quaternionToEulers = quaternionToEulers;
exports.append = append;
exports.floatToString = floatToString;
exports.getUniformTimes = getUniformTimes;
exports.clipTimesToClosestBVHTime = clipTimesToClosestBVHTime;
exports.distributeValue = distributeValue;
exports.lerpVector = lerpVector;
exports.lerpQuaternion = lerpQuaternion;
exports.lerpValues = lerpValues;
const quaternion_1 = __webpack_require__(/*! quaternion */ "./node_modules/quaternion/dist/quaternion.js");
exports.RAD_TO_DEG = 180.0 / Math.PI;
function toQuaternion(v) {
    const wSqr = 1.0 - (v.x * v.x + v.y * v.y + v.z * v.z);
    return new quaternion_1.Quaternion({ x: v.x, y: v.y, z: v.z, w: wSqr > 0 ? Math.sqrt(wSqr) : 0 });
}
function toEulers(truncatedQuanterion) {
    return toQuaternion({
        x: truncatedQuanterion.z,
        y: truncatedQuanterion.x,
        z: truncatedQuanterion.y
    }).toEuler().map(item1 => item1 * 180 / Math.PI);
}
function quaternionToEulers(quaternion) {
    const eulers = new quaternion_1.Quaternion({
        x: quaternion.z,
        y: quaternion.x,
        z: quaternion.y,
        w: quaternion.w
    }).toEuler().map(item1 => item1 * 180 / Math.PI);
    return { x: eulers[0], y: eulers[1], z: eulers[2] };
}
function append(text, times) {
    let result = "";
    for (let i = 0; i < times; i++) {
        result += text;
    }
    return result;
}
function floatToString(value, fraction) {
    return value.toLocaleString("un-US", { minimumFractionDigits: fraction }).replaceAll(",", "");
}
function lerpValue(x1, x2, t) {
    if (t > 1) {
        return x2;
    }
    if (t < 0) {
        return x1;
    }
    return x1 + (x2 - x1) * t;
}
function optimizedAmimationLength(duration, originalFrameLength) {
    if (originalFrameLength > duration) {
        return 2;
    }
    const hourLength = 3600 / originalFrameLength;
    let bestLength = 0;
    for (let i = 1; i < hourLength; i++) {
        const optimizedFrameLength = duration / i;
        const error = Math.abs(originalFrameLength - optimizedFrameLength);
        const prevError = Math.abs(originalFrameLength - bestLength);
        if (error < prevError) {
            bestLength = optimizedFrameLength;
        }
    }
    return bestLength;
}
function getUniformTimes(duration, singleFrameDuration) {
    const times = [];
    const optimizedFrameDuration = optimizedAmimationLength(duration, singleFrameDuration);
    const length = duration / optimizedFrameDuration;
    for (let i = 0; i < length; i++) {
        times.push(i * optimizedFrameDuration);
    }
    times[length - 1] = duration;
    return times;
}
function closest(left, right, value) {
    if (Math.abs(left - value) < Math.abs(right - value)) {
        return left;
    }
    return right;
}
function clipTimesToClosestBVHTime(animTimes, bvhTimes) {
    const fixedTimes = animTimes.map((item) => item);
    for (let i = 1; i < bvhTimes.length; i++) {
        for (let j = 0; j < animTimes.length; j++) {
            const animTime = animTimes[j];
            const bvhTimeLeft = bvhTimes[i - 1];
            const bvhTimeRight = bvhTimes[i];
            if ((bvhTimeLeft <= animTime) && (animTime <= bvhTimeRight)) {
                fixedTimes[j] = closest(bvhTimeLeft, bvhTimeRight, animTime);
            }
        }
    }
    return fixedTimes;
}
function getFactors(bvhTimes, animTimes) {
    return bvhTimes.map((item) => {
        if (item <= animTimes[0]) {
            return {
                leftAnimIndex: 0,
                rightAnimIndex: 1,
                factor: 0
            };
        }
        if (item >= animTimes[animTimes.length - 1]) {
            return {
                leftAnimIndex: animTimes.length - 2,
                rightAnimIndex: animTimes.length - 1,
                factor: 1
            };
        }
        for (let i = 1; i < animTimes.length; i++) {
            const leftTime = animTimes[i - 1];
            const rightTime = animTimes[i];
            if ((leftTime <= item) && (item < rightTime)) {
                const rangeSize = rightTime - leftTime;
                return {
                    leftAnimIndex: i - 1,
                    rightAnimIndex: i,
                    factor: (item - leftTime) / rangeSize
                };
            }
        }
    });
}
function distributeValue(position, rotation, channel, value) {
    const key = channel.toLowerCase()[0];
    const recipient = channel.includes("pos") ? position : rotation;
    recipient[key] = value;
}
function lerpVector(leftValue, rightValue, factor) {
    return {
        x: lerpValue(leftValue.x, rightValue.x, factor),
        y: lerpValue(leftValue.y, rightValue.y, factor),
        z: lerpValue(leftValue.z, rightValue.z, factor),
    };
}
function lerpQuaternion(leftValue, rightValue, factor) {
    return leftValue.slerp(rightValue)(factor);
}
function lerpValues(values, animTimes, uniformTimes, defaultValue, lerpFunction) {
    if (!(values === null || values === void 0 ? void 0 : values.length)) {
        return uniformTimes.map(item => defaultValue);
    }
    if (values.length == 1) {
        return uniformTimes.map(item => values[0]);
    }
    const factors = getFactors(uniformTimes, animTimes);
    return factors.map(item => {
        const leftFrame = values[item.leftAnimIndex];
        const rightFrame = values[item.rightAnimIndex];
        return lerpFunction(leftFrame, rightFrame, item.factor);
    });
}


/***/ })

/******/ 	});
/************************************************************************/
/******/ 	// The module cache
/******/ 	var __webpack_module_cache__ = {};
/******/ 	
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/ 		// Check if module is in cache
/******/ 		var cachedModule = __webpack_module_cache__[moduleId];
/******/ 		if (cachedModule !== undefined) {
/******/ 			return cachedModule.exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = __webpack_module_cache__[moduleId] = {
/******/ 			// no module.id needed
/******/ 			// no module.loaded needed
/******/ 			exports: {}
/******/ 		};
/******/ 	
/******/ 		// Execute the module function
/******/ 		__webpack_modules__[moduleId](module, module.exports, __webpack_require__);
/******/ 	
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/ 	
/************************************************************************/
var __webpack_exports__ = {};
// This entry needs to be wrapped in an IIFE because it needs to be isolated against other modules in the chunk.
(() => {
var exports = __webpack_exports__;
/*!**********************!*\
  !*** ./src/index.ts ***!
  \**********************/

Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.femaleOffsets = exports.maleOffsets = exports.visitNode = exports.serializeBVH = exports.toBVH = exports.parseBVH = exports.parseAnim = void 0;
const parse_1 = __webpack_require__(/*! ./parse */ "./src/parse.ts");
Object.defineProperty(exports, "parseAnim", ({ enumerable: true, get: function () { return parse_1.parseAnim; } }));
Object.defineProperty(exports, "parseBVH", ({ enumerable: true, get: function () { return parse_1.parseBVH; } }));
const convert_1 = __webpack_require__(/*! ./convert */ "./src/convert.ts");
Object.defineProperty(exports, "toBVH", ({ enumerable: true, get: function () { return convert_1.toBVH; } }));
Object.defineProperty(exports, "serializeBVH", ({ enumerable: true, get: function () { return convert_1.serializeBVH; } }));
Object.defineProperty(exports, "visitNode", ({ enumerable: true, get: function () { return convert_1.visitNode; } }));
const offsets_1 = __webpack_require__(/*! ./offsets */ "./src/offsets.ts");
Object.defineProperty(exports, "maleOffsets", ({ enumerable: true, get: function () { return offsets_1.maleOffsets; } }));
Object.defineProperty(exports, "femaleOffsets", ({ enumerable: true, get: function () { return offsets_1.femaleOffsets; } }));

})();

/******/ 	return __webpack_exports__;
/******/ })()
;
});
//# sourceMappingURL=data:application/json;charset=utf-8;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoiYnVuZGxlLmpzIiwibWFwcGluZ3MiOiJBQUFBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLENBQUM7QUFDRCxPOzs7Ozs7Ozs7O0FDVmE7O0FBRWI7QUFDQTtBQUNBO0FBQ0EsV0FBVyxRQUFRO0FBQ25CLFdBQVcsUUFBUTtBQUNuQixXQUFXLFFBQVE7QUFDbkIsV0FBVyxRQUFRO0FBQ25CO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0EsV0FBVyxRQUFRO0FBQ25CLFdBQVcsUUFBUTtBQUNuQixXQUFXLFFBQVE7QUFDbkIsV0FBVyxRQUFRO0FBQ25CO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQSxXQUFXLFFBQVE7QUFDbkIsV0FBVyxRQUFRO0FBQ25CLGFBQWE7QUFDYjtBQUNBOztBQUVBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTs7QUFFQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7O0FBRUE7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7O0FBRUEsbUJBQW1COztBQUVuQjtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQSxvQkFBb0IsbUJBQW1COztBQUV2QztBQUNBOztBQUVBO0FBQ0E7QUFDQSxRQUFRO0FBQ1I7QUFDQSxRQUFRO0FBQ1I7QUFDQSxRQUFROztBQUVSO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQSxZQUFZO0FBQ1o7QUFDQTs7QUFFQSxVQUFVOztBQUVWO0FBQ0E7QUFDQTs7QUFFQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsSUFBSTs7QUFFSjs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsTUFBTTtBQUNOO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTs7QUFFQTs7QUFFQTs7QUFFQTtBQUNBO0FBQ0EsTUFBTTtBQUNOO0FBQ0E7O0FBRUE7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQSxXQUFXLHNCQUFzQjtBQUNqQyxXQUFXLFNBQVM7QUFDcEIsV0FBVyxTQUFTO0FBQ3BCLFdBQVcsU0FBUztBQUNwQixhQUFhO0FBQ2I7QUFDQTs7QUFFQTtBQUNBO0FBQ0EsSUFBSTtBQUNKO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGFBQWEsc0JBQXNCO0FBQ25DLGFBQWEsU0FBUztBQUN0QixhQUFhLFNBQVM7QUFDdEIsYUFBYSxTQUFTO0FBQ3RCLGVBQWU7QUFDZjtBQUNBOztBQUVBOztBQUVBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsYUFBYSxzQkFBc0I7QUFDbkMsYUFBYSxTQUFTO0FBQ3RCLGFBQWEsU0FBUztBQUN0QixhQUFhLFNBQVM7QUFDdEIsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7O0FBRUE7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsR0FBRztBQUNIO0FBQ0E7QUFDQTtBQUNBLGVBQWU7QUFDZjtBQUNBOztBQUVBOztBQUVBO0FBQ0EsR0FBRztBQUNIO0FBQ0E7QUFDQTtBQUNBLGVBQWU7QUFDZjtBQUNBOztBQUVBOztBQUVBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0EsR0FBRztBQUNIO0FBQ0E7QUFDQTtBQUNBLGVBQWU7QUFDZjtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQTtBQUNBLGVBQWU7QUFDZjtBQUNBOztBQUVBOztBQUVBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTs7QUFFQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQTtBQUNBLGFBQWEsc0JBQXNCO0FBQ25DLGFBQWEsU0FBUztBQUN0QixhQUFhLFNBQVM7QUFDdEIsYUFBYSxTQUFTO0FBQ3RCLGVBQWU7QUFDZjtBQUNBOztBQUVBOztBQUVBOztBQUVBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsYUFBYSxRQUFRO0FBQ3JCLGVBQWU7QUFDZjtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsYUFBYSxzQkFBc0I7QUFDbkMsYUFBYSxTQUFTO0FBQ3RCLGFBQWEsU0FBUztBQUN0QixhQUFhLFNBQVM7QUFDdEIsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7O0FBRUE7O0FBRUE7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0E7QUFDQSxlQUFlO0FBQ2Y7QUFDQTs7QUFFQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7O0FBRUE7QUFDQSxpQ0FBaUM7QUFDakM7O0FBRUE7O0FBRUE7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsYUFBYSxzQkFBc0I7QUFDbkMsYUFBYSxTQUFTO0FBQ3RCLGFBQWEsU0FBUztBQUN0QixhQUFhLFNBQVM7QUFDdEIsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7O0FBRUE7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7O0FBRUE7QUFDQSxpQ0FBaUM7QUFDakM7O0FBRUE7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQSxlQUFlO0FBQ2Y7QUFDQTs7QUFFQTs7QUFFQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTs7QUFFQTtBQUNBOztBQUVBO0FBQ0EsR0FBRztBQUNIO0FBQ0E7QUFDQTtBQUNBLGVBQWU7QUFDZjtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsR0FBRztBQUNIO0FBQ0E7QUFDQTtBQUNBLGVBQWU7QUFDZjtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTs7QUFFQSxnREFBZ0Q7O0FBRWhEO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsYUFBYSxzQkFBc0I7QUFDbkMsYUFBYSxTQUFTO0FBQ3RCLGFBQWEsU0FBUztBQUN0QixhQUFhLFNBQVM7QUFDdEIsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7O0FBRUE7O0FBRUE7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTs7QUFFQTs7QUFFQTs7QUFFQTs7QUFFQSxZQUFZOztBQUVaO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsYUFBYSxzQkFBc0I7QUFDbkMsYUFBYSxTQUFTO0FBQ3RCLGFBQWEsU0FBUztBQUN0QixhQUFhLFNBQVM7QUFDdEIsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7O0FBRUE7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQSxlQUFlO0FBQ2Y7QUFDQTs7QUFFQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQSxlQUFlO0FBQ2Y7QUFDQTs7QUFFQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQSxlQUFlO0FBQ2Y7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7O0FBRUE7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsYUFBYSxVQUFVO0FBQ3ZCLGVBQWU7QUFDZjtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsYUFBYSxVQUFVO0FBQ3ZCLGVBQWU7QUFDZjtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsR0FBRztBQUNIO0FBQ0E7QUFDQTtBQUNBLGVBQWU7QUFDZjtBQUNBOztBQUVBOztBQUVBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0EsTUFBTTtBQUNOLGtDQUFrQztBQUNsQztBQUNBO0FBQ0EsR0FBRztBQUNIO0FBQ0E7QUFDQTtBQUNBLGVBQWU7QUFDZjtBQUNBOztBQUVBO0FBQ0EsNEJBQTRCOztBQUU1QiwwQkFBMEI7QUFDMUI7QUFDQSxxREFBcUQ7QUFDckQ7O0FBRUE7QUFDQSxvQ0FBb0M7QUFDcEM7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsYUFBYSxTQUFTO0FBQ3RCLGVBQWU7QUFDZjtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUEsK0NBQStDO0FBQy9DO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0Esa0NBQWtDO0FBQ2xDO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQSxtQ0FBbUM7QUFDbkM7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQSxlQUFlO0FBQ2Y7QUFDQTs7QUFFQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQTtBQUNBLGFBQWEsT0FBTztBQUNwQixlQUFlO0FBQ2Y7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLEdBQUc7O0FBRUg7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBOztBQUVBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBOztBQUVBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUEsZ0RBQWdEO0FBQ2hELCtDQUErQztBQUMvQztBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0EsV0FBVyxPQUFPO0FBQ2xCLFdBQVcsUUFBUTtBQUNuQixhQUFhO0FBQ2I7QUFDQTs7QUFFQTs7QUFFQTtBQUNBO0FBQ0E7O0FBRUE7O0FBRUE7QUFDQTs7QUFFQTs7QUFFQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsV0FBVyxPQUFPO0FBQ2xCLFdBQVcsT0FBTztBQUNsQjtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0EsTUFBTTtBQUNOO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQSxXQUFXLFFBQVE7QUFDbkIsV0FBVyxRQUFRO0FBQ25CLFdBQVcsUUFBUTtBQUNuQixXQUFXLFNBQVM7QUFDcEIsYUFBYTtBQUNiO0FBQ0E7O0FBRUE7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQSxXQUFXLFFBQVE7QUFDbkIsV0FBVyxRQUFRO0FBQ25CLFdBQVcsUUFBUTtBQUNuQixXQUFXLFNBQVM7QUFDcEIsYUFBYTtBQUNiO0FBQ0E7O0FBRUE7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQSw0Q0FBNEM7QUFDNUM7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUEseUJBQXlCO0FBQ3pCO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQSxXQUFXLE9BQU87QUFDbEIsYUFBYTtBQUNiO0FBQ0E7O0FBRUE7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQSxJQUFJO0FBQ0o7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQSw4QkFBOEI7O0FBRTlCOztBQUVBLGdCQUFnQjtBQUNoQjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxJQUFJO0FBQ0o7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsSUFBSTtBQUNKO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLElBQUk7QUFDSjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBLGtEQUFrRCxlQUFlO0FBQ2pFO0FBQ0E7QUFDQTs7Ozs7Ozs7Ozs7QUM5NENhO0FBQ2IsOENBQTZDLEVBQUUsYUFBYSxFQUFDO0FBQzdELGVBQWU7QUFDZixlQUFlO0FBQ2Y7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7Ozs7Ozs7Ozs7QUN6SWE7QUFDYiw4Q0FBNkMsRUFBRSxhQUFhLEVBQUM7QUFDN0QsaUJBQWlCO0FBQ2pCLG9CQUFvQjtBQUNwQixhQUFhO0FBQ2IsZ0JBQWdCLG1CQUFPLENBQUMsK0JBQVM7QUFDakMsb0JBQW9CLG1CQUFPLENBQUMsdUNBQWE7QUFDekMsa0JBQWtCLG1CQUFPLENBQUMsbUNBQVc7QUFDckMscUJBQXFCLG1CQUFPLENBQUMsZ0VBQVk7QUFDekM7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsMEVBQTBFO0FBQzFFO0FBQ0E7QUFDQSxzREFBc0Q7QUFDdEQ7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsMkNBQTJDLDBDQUEwQztBQUNyRjtBQUNBLHVCQUF1QjtBQUN2QjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsMkJBQTJCLGdCQUFnQjtBQUMzQztBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSw4Q0FBOEMsWUFBWSxrS0FBa0s7QUFDNU47QUFDQTtBQUNBO0FBQ0E7QUFDQSw4Q0FBOEMsWUFBWSxrS0FBa0s7QUFDNU47QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsYUFBYTtBQUNiO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0Esd0JBQXdCO0FBQ3hCO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpSkFBaUosa0JBQWtCO0FBQ25LO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxhQUFhO0FBQ2I7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0Esb0JBQW9CLDZCQUE2QjtBQUNqRDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7Ozs7Ozs7Ozs7O0FDdExhO0FBQ2IsOENBQTZDLEVBQUUsYUFBYSxFQUFDO0FBQzdELGlCQUFpQjtBQUNqQixpQkFBaUI7QUFDakI7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EseUVBQXlFO0FBQ3pFO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EseUVBQXlFO0FBQ3pFO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EseUVBQXlFO0FBQ3pFO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUZBQWlGO0FBQ2pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUZBQWlGO0FBQ2pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUZBQWlGO0FBQ2pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUZBQWlGO0FBQ2pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUZBQWlGO0FBQ2pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUZBQWlGO0FBQ2pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUZBQWlGO0FBQ2pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUZBQWlGO0FBQ2pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUZBQWlGO0FBQ2pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUZBQWlGO0FBQ2pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUZBQWlGO0FBQ2pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUZBQWlGO0FBQ2pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUZBQWlGO0FBQ2pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUZBQWlGO0FBQ2pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlGQUFpRjtBQUNqRjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSx5RkFBeUY7QUFDekY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpR0FBaUc7QUFDakc7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpR0FBaUc7QUFDakc7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpR0FBaUc7QUFDakc7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSx5RkFBeUY7QUFDekY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSx5RkFBeUY7QUFDekY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSx5RkFBeUY7QUFDekY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSx5RkFBeUY7QUFDekY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUZBQWlGO0FBQ2pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUZBQWlGO0FBQ2pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUZBQWlGO0FBQ2pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EseURBQXlEO0FBQ3pEO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSx5RkFBeUY7QUFDekY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHlGQUF5RjtBQUN6RjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EseUZBQXlGO0FBQ3pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSx5RkFBeUY7QUFDekY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSx5REFBeUQ7QUFDekQ7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHlGQUF5RjtBQUN6RjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EseUZBQXlGO0FBQ3pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSx5RkFBeUY7QUFDekY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHlGQUF5RjtBQUN6RjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHlEQUF5RDtBQUN6RDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHlGQUF5RjtBQUN6RjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlFQUFpRTtBQUNqRTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHlGQUF5RjtBQUN6RjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxTQUFTO0FBQ1Q7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsU0FBUztBQUNUO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLFNBQVM7QUFDVDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsU0FBUztBQUNUO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsU0FBUztBQUNUO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUJBQWlCO0FBQ2pCO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7Ozs7Ozs7Ozs7QUNyMUJhO0FBQ2IsOENBQTZDLEVBQUUsYUFBYSxFQUFDO0FBQzdELG1CQUFtQixHQUFHLHFCQUFxQjtBQUMzQyxxQkFBcUI7QUFDckI7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsbUJBQW1CO0FBQ25CO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7Ozs7Ozs7Ozs7QUNwNkRhO0FBQ2IsOENBQTZDLEVBQUUsYUFBYSxFQUFDO0FBQzdELGlCQUFpQjtBQUNqQiw2QkFBNkI7QUFDN0IsZ0JBQWdCO0FBQ2hCLGdCQUFnQixtQkFBTyxDQUFDLCtCQUFTO0FBQ2pDO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsdUJBQXVCO0FBQ3ZCO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHdCQUF3QixVQUFVO0FBQ2xDO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0Esb0JBQW9CLGdCQUFnQjtBQUNwQztBQUNBO0FBQ0E7QUFDQTtBQUNBLHdCQUF3QixrQkFBa0I7QUFDMUM7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpQ0FBaUMsNEhBQTRIO0FBQzdKO0FBQ0E7QUFDQTtBQUNBLHdCQUF3QixrQkFBa0I7QUFDMUM7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpQ0FBaUMsNEhBQTRIO0FBQzdKO0FBQ0Esc0JBQXNCLDBEQUEwRDtBQUNoRjtBQUNBO0FBQ0E7QUFDQSxvQkFBb0IscUJBQXFCO0FBQ3pDO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLFNBQVM7QUFDVDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0wsYUFBYTtBQUNiO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMLHVDQUF1QywrRUFBK0U7QUFDdEgsdUNBQXVDLCtFQUErRTtBQUN0SDtBQUNBO0FBQ0E7QUFDQSwwQkFBMEI7QUFDMUI7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxvQkFBb0IsdUJBQXVCO0FBQzNDO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlCQUFpQjtBQUNqQjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsYUFBYTtBQUNiO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsdUJBQXVCO0FBQ3ZCLHVCQUF1QjtBQUN2QjtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0EsK0JBQStCLG9CQUFvQjtBQUNuRDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlDQUFpQyxrQkFBa0I7QUFDbkQsNEJBQTRCO0FBQzVCO0FBQ0EsK0NBQStDO0FBQy9DO0FBQ0E7QUFDQSxLQUFLLFNBQVM7QUFDZCw0QkFBNEI7QUFDNUI7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLLFNBQVM7QUFDZDtBQUNBO0FBQ0E7QUFDQSwyQ0FBMkMsZUFBZTtBQUMxRDtBQUNBO0FBQ0Esb0JBQW9CLDRCQUE0QjtBQUNoRDtBQUNBO0FBQ0E7QUFDQTtBQUNBOzs7Ozs7Ozs7OztBQ2pQYTtBQUNiLDhDQUE2QyxFQUFFLGFBQWEsRUFBQztBQUM3RCxrQkFBa0I7QUFDbEIsb0JBQW9CO0FBQ3BCLGdCQUFnQjtBQUNoQiwwQkFBMEI7QUFDMUIsY0FBYztBQUNkLHFCQUFxQjtBQUNyQix1QkFBdUI7QUFDdkIsaUNBQWlDO0FBQ2pDLHVCQUF1QjtBQUN2QixrQkFBa0I7QUFDbEIsc0JBQXNCO0FBQ3RCLGtCQUFrQjtBQUNsQixxQkFBcUIsbUJBQU8sQ0FBQyxnRUFBWTtBQUN6QyxrQkFBa0I7QUFDbEI7QUFDQTtBQUNBLHlDQUF5QywyREFBMkQ7QUFDcEc7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMLGFBQWE7QUFDYjtBQUNBO0FBQ0E7QUFDQSxvQkFBb0IsV0FBVztBQUMvQjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsMkNBQTJDLGlDQUFpQztBQUM1RTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLG9CQUFvQixnQkFBZ0I7QUFDcEM7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxvQkFBb0IsWUFBWTtBQUNoQztBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLG9CQUFvQixxQkFBcUI7QUFDekMsd0JBQXdCLHNCQUFzQjtBQUM5QztBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0Esd0JBQXdCLHNCQUFzQjtBQUM5QztBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMOzs7Ozs7O1VDL0pBO1VBQ0E7O1VBRUE7VUFDQTtVQUNBO1VBQ0E7VUFDQTtVQUNBO1VBQ0E7VUFDQTtVQUNBO1VBQ0E7VUFDQTtVQUNBO1VBQ0E7O1VBRUE7VUFDQTs7VUFFQTtVQUNBO1VBQ0E7Ozs7Ozs7Ozs7QUN0QmE7QUFDYiw4Q0FBNkMsRUFBRSxhQUFhLEVBQUM7QUFDN0QscUJBQXFCLEdBQUcsbUJBQW1CLEdBQUcsaUJBQWlCLEdBQUcsb0JBQW9CLEdBQUcsYUFBYSxHQUFHLGdCQUFnQixHQUFHLGlCQUFpQjtBQUM3SSxnQkFBZ0IsbUJBQU8sQ0FBQywrQkFBUztBQUNqQyw2Q0FBNEMsRUFBRSxxQ0FBcUMsNkJBQTZCLEVBQUM7QUFDakgsNENBQTJDLEVBQUUscUNBQXFDLDRCQUE0QixFQUFDO0FBQy9HLGtCQUFrQixtQkFBTyxDQUFDLG1DQUFXO0FBQ3JDLHlDQUF3QyxFQUFFLHFDQUFxQywyQkFBMkIsRUFBQztBQUMzRyxnREFBK0MsRUFBRSxxQ0FBcUMsa0NBQWtDLEVBQUM7QUFDekgsNkNBQTRDLEVBQUUscUNBQXFDLCtCQUErQixFQUFDO0FBQ25ILGtCQUFrQixtQkFBTyxDQUFDLG1DQUFXO0FBQ3JDLCtDQUE4QyxFQUFFLHFDQUFxQyxpQ0FBaUMsRUFBQztBQUN2SCxpREFBZ0QsRUFBRSxxQ0FBcUMsbUNBQW1DLEVBQUMiLCJzb3VyY2VzIjpbIndlYnBhY2s6Ly9BbmltVG9Cdmgvd2VicGFjay91bml2ZXJzYWxNb2R1bGVEZWZpbml0aW9uIiwid2VicGFjazovL0FuaW1Ub0J2aC8uL25vZGVfbW9kdWxlcy9xdWF0ZXJuaW9uL2Rpc3QvcXVhdGVybmlvbi5qcyIsIndlYnBhY2s6Ly9BbmltVG9CdmgvLi9zcmMvYWxpYXNlcy50cyIsIndlYnBhY2s6Ly9BbmltVG9CdmgvLi9zcmMvY29udmVydC50cyIsIndlYnBhY2s6Ly9BbmltVG9CdmgvLi9zcmMvaGllcmFyY2h5LnRzIiwid2VicGFjazovL0FuaW1Ub0J2aC8uL3NyYy9vZmZzZXRzLnRzIiwid2VicGFjazovL0FuaW1Ub0J2aC8uL3NyYy9wYXJzZS50cyIsIndlYnBhY2s6Ly9BbmltVG9CdmgvLi9zcmMvdXRpbHMudHMiLCJ3ZWJwYWNrOi8vQW5pbVRvQnZoL3dlYnBhY2svYm9vdHN0cmFwIiwid2VicGFjazovL0FuaW1Ub0J2aC8uL3NyYy9pbmRleC50cyJdLCJzb3VyY2VzQ29udGVudCI6WyIoZnVuY3Rpb24gd2VicGFja1VuaXZlcnNhbE1vZHVsZURlZmluaXRpb24ocm9vdCwgZmFjdG9yeSkge1xuXHRpZih0eXBlb2YgZXhwb3J0cyA9PT0gJ29iamVjdCcgJiYgdHlwZW9mIG1vZHVsZSA9PT0gJ29iamVjdCcpXG5cdFx0bW9kdWxlLmV4cG9ydHMgPSBmYWN0b3J5KCk7XG5cdGVsc2UgaWYodHlwZW9mIGRlZmluZSA9PT0gJ2Z1bmN0aW9uJyAmJiBkZWZpbmUuYW1kKVxuXHRcdGRlZmluZShbXSwgZmFjdG9yeSk7XG5cdGVsc2UgaWYodHlwZW9mIGV4cG9ydHMgPT09ICdvYmplY3QnKVxuXHRcdGV4cG9ydHNbXCJBbmltVG9CdmhcIl0gPSBmYWN0b3J5KCk7XG5cdGVsc2Vcblx0XHRyb290W1wiQW5pbVRvQnZoXCJdID0gZmFjdG9yeSgpO1xufSkodHlwZW9mIHNlbGYgIT09ICd1bmRlZmluZWQnID8gc2VsZiA6IHRoaXMsICgpID0+IHtcbnJldHVybiAiLCIndXNlIHN0cmljdCc7XG5cbi8qKlxuICogQ3JlYXRlcyBhIG5ldyBRdWF0ZXJuaW9uIG9iamVjdFxuICogXG4gKiBAcGFyYW0ge251bWJlcn0gdyBcbiAqIEBwYXJhbSB7bnVtYmVyfSB4IFxuICogQHBhcmFtIHtudW1iZXJ9IHkgXG4gKiBAcGFyYW0ge251bWJlcn0geiBcbiAqIEByZXR1cm5zIFxuICovXG5mdW5jdGlvbiBuZXdRdWF0ZXJuaW9uKHcsIHgsIHksIHopIHtcbiAgY29uc3QgZiA9IE9iamVjdC5jcmVhdGUoUXVhdGVybmlvbi5wcm90b3R5cGUpO1xuXG4gIGZbJ3cnXSA9IHc7XG4gIGZbJ3gnXSA9IHg7XG4gIGZbJ3knXSA9IHk7XG4gIGZbJ3onXSA9IHo7XG5cbiAgcmV0dXJuIGY7XG59XG5cbi8qKlxuICogQ3JlYXRlcyBhIG5ldyBub3JtYWxpemVkIFF1YXRlcm5pb24gb2JqZWN0XG4gKlxuICogQHBhcmFtIHtudW1iZXJ9IHdcbiAqIEBwYXJhbSB7bnVtYmVyfSB4XG4gKiBAcGFyYW0ge251bWJlcn0geVxuICogQHBhcmFtIHtudW1iZXJ9IHpcbiAqIEByZXR1cm5zXG4gKi9cbmZ1bmN0aW9uIG5ld05vcm1hbGl6ZWQodywgeCwgeSwgeikge1xuICBjb25zdCBmID0gT2JqZWN0LmNyZWF0ZShRdWF0ZXJuaW9uLnByb3RvdHlwZSk7XG5cbiAgLy8gV2UgYXNzdW1lIHxRfCA+IDAgZm9yIGludGVybmFsIHVzYWdlXG4gIGNvbnN0IGlsID0gMSAvIE1hdGguc3FydCh3ICogdyArIHggKiB4ICsgeSAqIHkgKyB6ICogeik7XG5cbiAgZlsndyddID0gdyAqIGlsO1xuICBmWyd4J10gPSB4ICogaWw7XG4gIGZbJ3knXSA9IHkgKiBpbDtcbiAgZlsneiddID0geiAqIGlsO1xuXG4gIHJldHVybiBmO1xufVxuXG4vKipcbiAqIENhbGN1bGF0ZXMgbG9nKHNxcnQoYV4yK2JeMikpIGluIGEgd2F5IHRvIGF2b2lkIG92ZXJmbG93c1xuICpcbiAqIEBwYXJhbSB7bnVtYmVyfSBhXG4gKiBAcGFyYW0ge251bWJlcn0gYlxuICogQHJldHVybnMge251bWJlcn1cbiAqL1xuZnVuY3Rpb24gbG9nSHlwb3QoYSwgYikge1xuXG4gIGNvbnN0IF9hID0gTWF0aC5hYnMoYSk7XG4gIGNvbnN0IF9iID0gTWF0aC5hYnMoYik7XG5cbiAgaWYgKGEgPT09IDApIHtcbiAgICByZXR1cm4gTWF0aC5sb2coX2IpO1xuICB9XG5cbiAgaWYgKGIgPT09IDApIHtcbiAgICByZXR1cm4gTWF0aC5sb2coX2EpO1xuICB9XG5cbiAgaWYgKF9hIDwgMzAwMCAmJiBfYiA8IDMwMDApIHtcbiAgICByZXR1cm4gMC41ICogTWF0aC5sb2coYSAqIGEgKyBiICogYik7XG4gIH1cblxuICBhID0gYSAvIDI7XG4gIGIgPSBiIC8gMjtcblxuICByZXR1cm4gMC41ICogTWF0aC5sb2coYSAqIGEgKyBiICogYikgKyBNYXRoLkxOMjtcbn1cblxuLypcbiAqIFRlbXBvcmFyeSBwYXJzaW5nIG9iamVjdCB0byBhdm9pZCByZS1hbGxvY2F0aW9uc1xuICpcbiAqL1xuY29uc3QgUCA9IE9iamVjdC5jcmVhdGUoUXVhdGVybmlvbi5wcm90b3R5cGUpO1xuXG5mdW5jdGlvbiBwYXJzZShkZXN0LCB3LCB4LCB5LCB6KSB7XG5cbiAgLy8gTW9zdCBjb21tb24gaW50ZXJuYWwgdXNlIGNhc2Ugd2l0aCA0IHBhcmFtc1xuICBpZiAoeiAhPT0gdW5kZWZpbmVkKSB7XG4gICAgZGVzdFsndyddID0gdztcbiAgICBkZXN0Wyd4J10gPSB4O1xuICAgIGRlc3RbJ3knXSA9IHk7XG4gICAgZGVzdFsneiddID0gejtcbiAgICByZXR1cm47XG4gIH1cblxuICBpZiAodHlwZW9mIHcgPT09ICdvYmplY3QnICYmIHkgPT09IHVuZGVmaW5lZCkge1xuXG4gICAgLy8gQ2hlY2sgZm9yIHF1YXRzLCBmb3IgZXhhbXBsZSB3aGVuIGFuIG9iamVjdCBnZXRzIGNsb25lZFxuICAgIGlmICgndycgaW4gdyB8fCAneCcgaW4gdyB8fCAneScgaW4gdyB8fCAneicgaW4gdykge1xuICAgICAgZGVzdFsndyddID0gd1sndyddIHx8IDA7XG4gICAgICBkZXN0Wyd4J10gPSB3Wyd4J10gfHwgMDtcbiAgICAgIGRlc3RbJ3knXSA9IHdbJ3knXSB8fCAwO1xuICAgICAgZGVzdFsneiddID0gd1sneiddIHx8IDA7XG4gICAgICByZXR1cm47XG4gICAgfVxuXG4gICAgLy8gQ2hlY2sgZm9yIGNvbXBsZXggbnVtYmVyc1xuICAgIGlmICgncmUnIGluIHcgJiYgJ2ltJyBpbiB3KSB7XG4gICAgICBkZXN0Wyd3J10gPSB3WydyZSddO1xuICAgICAgZGVzdFsneCddID0gd1snaW0nXTtcbiAgICAgIGRlc3RbJ3knXSA9IDA7XG4gICAgICBkZXN0Wyd6J10gPSAwO1xuICAgICAgcmV0dXJuO1xuICAgIH1cblxuICAgIC8vIENoZWNrIGZvciBhcnJheVxuICAgIGlmICh3Lmxlbmd0aCA9PT0gNCkge1xuICAgICAgZGVzdFsndyddID0gd1swXTtcbiAgICAgIGRlc3RbJ3gnXSA9IHdbMV07XG4gICAgICBkZXN0Wyd5J10gPSB3WzJdO1xuICAgICAgZGVzdFsneiddID0gd1szXTtcbiAgICAgIHJldHVybjtcbiAgICB9XG5cbiAgICAvLyBDaGVjayBmb3IgYXVnbWVudGVkIHZlY3RvclxuICAgIGlmICh3Lmxlbmd0aCA9PT0gMykge1xuICAgICAgZGVzdFsndyddID0gMDtcbiAgICAgIGRlc3RbJ3gnXSA9IHdbMF07XG4gICAgICBkZXN0Wyd5J10gPSB3WzFdO1xuICAgICAgZGVzdFsneiddID0gd1syXTtcbiAgICAgIHJldHVybjtcbiAgICB9XG5cbiAgICB0aHJvdyBuZXcgRXJyb3IoJ0ludmFsaWQgb2JqZWN0Jyk7XG4gIH1cblxuICAvLyBQYXJzZSBzdHJpbmcgdmFsdWVzXG4gIGlmICh0eXBlb2YgdyA9PT0gJ3N0cmluZycgJiYgeSA9PT0gdW5kZWZpbmVkKSB7XG5cbiAgICBjb25zdCB0b2tlbnMgPSB3Lm1hdGNoKC9cXGQrXFwuP1xcZCplWystXT9cXGQrfFxcZCtcXC4/XFxkKnxcXC5cXGQrfC4vZyk7XG4gICAgbGV0IHBsdXMgPSAxO1xuICAgIGxldCBtaW51cyA9IDA7XG5cbiAgICBjb25zdCBpTWFwID0geyAnaSc6ICd4JywgJ2onOiAneScsICdrJzogJ3onIH07XG5cbiAgICBpZiAodG9rZW5zID09PSBudWxsKSB7XG4gICAgICB0aHJvdyBuZXcgRXJyb3IoJ1BhcnNlIGVycm9yJyk7XG4gICAgfVxuXG4gICAgLy8gUmVzZXQgdGhlIGN1cnJlbnQgc3RhdGVcbiAgICBkZXN0Wyd3J10gPVxuICAgICAgZGVzdFsneCddID1cbiAgICAgIGRlc3RbJ3knXSA9XG4gICAgICBkZXN0Wyd6J10gPSAwO1xuXG4gICAgZm9yIChsZXQgaSA9IDA7IGkgPCB0b2tlbnMubGVuZ3RoOyBpKyspIHtcblxuICAgICAgbGV0IGMgPSB0b2tlbnNbaV07XG4gICAgICBsZXQgZCA9IHRva2Vuc1tpICsgMV07XG5cbiAgICAgIGlmIChjID09PSAnICcgfHwgYyA9PT0gJ1xcdCcgfHwgYyA9PT0gJ1xcbicpIHtcbiAgICAgICAgLyogdm9pZCAqL1xuICAgICAgfSBlbHNlIGlmIChjID09PSAnKycpIHtcbiAgICAgICAgcGx1cysrO1xuICAgICAgfSBlbHNlIGlmIChjID09PSAnLScpIHtcbiAgICAgICAgbWludXMrKztcbiAgICAgIH0gZWxzZSB7XG5cbiAgICAgICAgaWYgKHBsdXMgKyBtaW51cyA9PT0gMCkge1xuICAgICAgICAgIHRocm93IG5ldyBFcnJvcignUGFyc2UgZXJyb3InICsgYyk7XG4gICAgICAgIH1cbiAgICAgICAgbGV0IGcgPSBpTWFwW2NdO1xuXG4gICAgICAgIC8vIElzIHRoZSBjdXJyZW50IHRva2VuIGFuIGltYWdpbmFyeSBzaWduP1xuICAgICAgICBpZiAoZyAhPT0gdW5kZWZpbmVkKSB7XG5cbiAgICAgICAgICAvLyBJcyB0aGUgZm9sbG93aW5nIHRva2VuIGEgbnVtYmVyP1xuICAgICAgICAgIGlmIChkICE9PSAnICcgJiYgIWlzTmFOKGQpKSB7XG4gICAgICAgICAgICBjID0gZDtcbiAgICAgICAgICAgIGkrKztcbiAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgYyA9ICcxJztcbiAgICAgICAgICB9XG5cbiAgICAgICAgfSBlbHNlIHtcblxuICAgICAgICAgIGlmIChpc05hTihjKSkge1xuICAgICAgICAgICAgdGhyb3cgbmV3IEVycm9yKCdQYXJzZXIgZXJyb3InKTtcbiAgICAgICAgICB9XG5cbiAgICAgICAgICBnID0gaU1hcFtkXTtcblxuICAgICAgICAgIGlmIChnICE9PSB1bmRlZmluZWQpIHtcbiAgICAgICAgICAgIGkrKztcbiAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICBkZXN0W2cgfHwgJ3cnXSArPSBwYXJzZUZsb2F0KChtaW51cyAlIDIgPyAnLScgOiAnJykgKyBjKTtcbiAgICAgICAgcGx1cyA9IG1pbnVzID0gMDtcbiAgICAgIH1cbiAgICB9XG5cbiAgICAvLyBTdGlsbCBzb21ldGhpbmcgb24gdGhlIHN0YWNrXG4gICAgaWYgKHBsdXMgKyBtaW51cyA+IDApIHtcbiAgICAgIHRocm93IG5ldyBFcnJvcignUGFyc2VyIGVycm9yJyk7XG4gICAgfVxuICAgIHJldHVybjtcbiAgfVxuXG4gIC8vIElmIG5vIHNpbmdsZSB2YXJpYWJsZSB3YXMgZ2l2ZW4gQU5EIGl0IHdhcyB0aGUgY29uc3RydWN0b3IsIHNldCBpdCB0byB0aGUgaWRlbnRpdHlcbiAgaWYgKHcgPT09IHVuZGVmaW5lZCAmJiBkZXN0ICE9PSBQKSB7XG4gICAgZGVzdFsndyddID0gMTtcbiAgICBkZXN0Wyd4J10gPVxuICAgICAgZGVzdFsneSddID1cbiAgICAgIGRlc3RbJ3onXSA9IDA7XG4gIH0gZWxzZSB7XG5cbiAgICBkZXN0Wyd3J10gPSB3IHx8IDA7XG5cbiAgICAvLyBOb3RlOiBUaGlzIGlzbid0IGZyb21BeGlzKCksIGl0J3MganVzdCBzeW50YWN0aWMgc3VnYXIhXG4gICAgaWYgKHggJiYgeC5sZW5ndGggPT09IDMpIHtcbiAgICAgIGRlc3RbJ3gnXSA9IHhbMF07XG4gICAgICBkZXN0Wyd5J10gPSB4WzFdO1xuICAgICAgZGVzdFsneiddID0geFsyXTtcbiAgICB9IGVsc2Uge1xuICAgICAgZGVzdFsneCddID0geCB8fCAwO1xuICAgICAgZGVzdFsneSddID0geSB8fCAwO1xuICAgICAgZGVzdFsneiddID0geiB8fCAwO1xuICAgIH1cbiAgfVxufVxuXG5mdW5jdGlvbiBudW1Ub1N0cihuLCBjaGFyLCBwcmV2KSB7XG5cbiAgbGV0IHJldCA9ICcnO1xuXG4gIGlmIChuICE9PSAwKSB7XG5cbiAgICBpZiAocHJldiAhPT0gJycpIHtcbiAgICAgIHJldCArPSBuIDwgMCA/ICcgLSAnIDogJyArICc7XG4gICAgfSBlbHNlIGlmIChuIDwgMCkge1xuICAgICAgcmV0ICs9ICctJztcbiAgICB9XG5cbiAgICBuID0gTWF0aC5hYnMobik7XG5cbiAgICBpZiAoMSAhPT0gbiB8fCBjaGFyID09PSAnJykge1xuICAgICAgcmV0ICs9IG47XG4gICAgfVxuICAgIHJldCArPSBjaGFyO1xuICB9XG4gIHJldHVybiByZXQ7XG59XG5cbi8qKlxuICogUXVhdGVybmlvbiBjb25zdHJ1Y3RvclxuICpcbiAqIEBjb25zdHJ1Y3RvclxuICogQHBhcmFtIHtudW1iZXJ8T2JqZWN0fHN0cmluZ30gdyByZWFsXG4gKiBAcGFyYW0ge251bWJlcj19IHggaW1hZ1xuICogQHBhcmFtIHtudW1iZXI9fSB5IGltYWdcbiAqIEBwYXJhbSB7bnVtYmVyPX0geiBpbWFnXG4gKiBAcmV0dXJucyB7UXVhdGVybmlvbn1cbiAqL1xuZnVuY3Rpb24gUXVhdGVybmlvbih3LCB4LCB5LCB6KSB7XG5cbiAgaWYgKHRoaXMgaW5zdGFuY2VvZiBRdWF0ZXJuaW9uKSB7XG4gICAgcGFyc2UodGhpcywgdywgeCwgeSwgeik7XG4gIH0gZWxzZSB7XG4gICAgY29uc3QgdCA9IE9iamVjdC5jcmVhdGUoUXVhdGVybmlvbi5wcm90b3R5cGUpO1xuICAgIHBhcnNlKHQsIHcsIHgsIHksIHopO1xuICAgIHJldHVybiB0O1xuICB9XG59XG5cblF1YXRlcm5pb24ucHJvdG90eXBlID0ge1xuICAndyc6IDEsXG4gICd4JzogMCxcbiAgJ3knOiAwLFxuICAneic6IDAsXG4gIC8qKlxuICAgKiBBZGRzIHR3byBxdWF0ZXJuaW9ucyBRMSBhbmQgUTJcbiAgICpcbiAgICogQHBhcmFtIHtudW1iZXJ8T2JqZWN0fHN0cmluZ30gdyByZWFsXG4gICAqIEBwYXJhbSB7bnVtYmVyPX0geCBpbWFnXG4gICAqIEBwYXJhbSB7bnVtYmVyPX0geSBpbWFnXG4gICAqIEBwYXJhbSB7bnVtYmVyPX0geiBpbWFnXG4gICAqIEByZXR1cm5zIHtRdWF0ZXJuaW9ufVxuICAgKi9cbiAgJ2FkZCc6IGZ1bmN0aW9uICh3LCB4LCB5LCB6KSB7XG5cbiAgICBwYXJzZShQLCB3LCB4LCB5LCB6KTtcblxuICAgIC8vIFExICsgUTIgOj0gW3cxLCB2MV0gKyBbdzIsIHYyXSA9IFt3MSArIHcyLCB2MSArIHYyXVxuXG4gICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oXG4gICAgICB0aGlzWyd3J10gKyBQWyd3J10sXG4gICAgICB0aGlzWyd4J10gKyBQWyd4J10sXG4gICAgICB0aGlzWyd5J10gKyBQWyd5J10sXG4gICAgICB0aGlzWyd6J10gKyBQWyd6J10pO1xuICB9LFxuICAvKipcbiAgICogU3VidHJhY3RzIGEgcXVhdGVybmlvbnMgUTIgZnJvbSBRMVxuICAgKlxuICAgKiBAcGFyYW0ge251bWJlcnxPYmplY3R8c3RyaW5nfSB3IHJlYWxcbiAgICogQHBhcmFtIHtudW1iZXI9fSB4IGltYWdcbiAgICogQHBhcmFtIHtudW1iZXI9fSB5IGltYWdcbiAgICogQHBhcmFtIHtudW1iZXI9fSB6IGltYWdcbiAgICogQHJldHVybnMge1F1YXRlcm5pb259XG4gICAqL1xuICAnc3ViJzogZnVuY3Rpb24gKHcsIHgsIHksIHopIHtcblxuICAgIHBhcnNlKFAsIHcsIHgsIHksIHopO1xuXG4gICAgLy8gUTEgLSBRMiA6PSBRMSArICgtUTIpXG4gICAgLy8gICAgICAgICAgPSBbdzEsIHYxXSAtIFt3MiwgdjJdID0gW3cxIC0gdzIsIHYxIC0gdjJdXG5cbiAgICByZXR1cm4gbmV3UXVhdGVybmlvbihcbiAgICAgIHRoaXNbJ3cnXSAtIFBbJ3cnXSxcbiAgICAgIHRoaXNbJ3gnXSAtIFBbJ3gnXSxcbiAgICAgIHRoaXNbJ3knXSAtIFBbJ3knXSxcbiAgICAgIHRoaXNbJ3onXSAtIFBbJ3onXSk7XG4gIH0sXG4gIC8qKlxuICAgKiBDYWxjdWxhdGVzIHRoZSBhZGRpdGl2ZSBpbnZlcnNlLCBvciBzaW1wbHkgaXQgbmVnYXRlcyB0aGUgcXVhdGVybmlvblxuICAgKlxuICAgKiBAcmV0dXJucyB7UXVhdGVybmlvbn1cbiAgICovXG4gICduZWcnOiBmdW5jdGlvbiAoKSB7XG5cbiAgICAvLyAtUSA6PSBbLXcsIC12XVxuXG4gICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oLXRoaXNbJ3cnXSwgLXRoaXNbJ3gnXSwgLXRoaXNbJ3knXSwgLXRoaXNbJ3onXSk7XG4gIH0sXG4gIC8qKlxuICAgKiBDYWxjdWxhdGVzIHRoZSBsZW5ndGgvbW9kdWx1cy9tYWduaXR1ZGUgb3IgdGhlIG5vcm0gb2YgYSBxdWF0ZXJuaW9uXG4gICAqXG4gICAqIEByZXR1cm5zIHtudW1iZXJ9XG4gICAqL1xuICAnbm9ybSc6IGZ1bmN0aW9uICgpIHtcblxuICAgIC8vIHxRfCA6PSBzcXJ0KHxRfF4yKVxuXG4gICAgLy8gVGhlIHVuaXQgcXVhdGVybmlvbiBoYXMgfFF8ID0gMVxuXG4gICAgY29uc3QgdyA9IHRoaXNbJ3cnXTtcbiAgICBjb25zdCB4ID0gdGhpc1sneCddO1xuICAgIGNvbnN0IHkgPSB0aGlzWyd5J107XG4gICAgY29uc3QgeiA9IHRoaXNbJ3onXTtcblxuICAgIHJldHVybiBNYXRoLnNxcnQodyAqIHcgKyB4ICogeCArIHkgKiB5ICsgeiAqIHopO1xuICB9LFxuICAvKipcbiAgICogQ2FsY3VsYXRlcyB0aGUgc3F1YXJlZCBsZW5ndGgvbW9kdWx1cy9tYWduaXR1ZGUgb3IgdGhlIG5vcm0gb2YgYSBxdWF0ZXJuaW9uXG4gICAqXG4gICAqIEByZXR1cm5zIHtudW1iZXJ9XG4gICAqL1xuICAnbm9ybVNxJzogZnVuY3Rpb24gKCkge1xuXG4gICAgLy8gfFF8XjIgOj0gW3csIHZdICogW3csIC12XVxuICAgIC8vICAgICAgICA9IFt3XjIgKyBkb3QodiwgdiksIC13ICogdiArIHcgKiB2ICsgY3Jvc3ModiwgLXYpXVxuICAgIC8vICAgICAgICA9IFt3XjIgKyB8dnxeMiwgMF1cbiAgICAvLyAgICAgICAgPSBbd14yICsgZG90KHYsIHYpLCAwXVxuICAgIC8vICAgICAgICA9IGRvdChRLCBRKVxuICAgIC8vICAgICAgICA9IFEgKiBRJ1xuXG4gICAgY29uc3QgdyA9IHRoaXNbJ3cnXTtcbiAgICBjb25zdCB4ID0gdGhpc1sneCddO1xuICAgIGNvbnN0IHkgPSB0aGlzWyd5J107XG4gICAgY29uc3QgeiA9IHRoaXNbJ3onXTtcblxuICAgIHJldHVybiB3ICogdyArIHggKiB4ICsgeSAqIHkgKyB6ICogejtcbiAgfSxcbiAgLyoqXG4gICAqIE5vcm1hbGl6ZXMgdGhlIHF1YXRlcm5pb24gdG8gaGF2ZSB8UXwgPSAxIGFzIGxvbmcgYXMgdGhlIG5vcm0gaXMgbm90IHplcm9cbiAgICogQWx0ZXJuYXRpdmUgbmFtZXMgYXJlIHRoZSBzaWdudW0sIHVuaXQgb3IgdmVyc29yXG4gICAqXG4gICAqIEByZXR1cm5zIHtRdWF0ZXJuaW9ufVxuICAgKi9cbiAgJ25vcm1hbGl6ZSc6IGZ1bmN0aW9uICgpIHtcblxuICAgIC8vIFEqIDo9IFEgLyB8UXxcblxuICAgIC8vIHVucm9sbGVkIFEuc2NhbGUoMSAvIFEubm9ybSgpKVxuXG4gICAgY29uc3QgdyA9IHRoaXNbJ3cnXTtcbiAgICBjb25zdCB4ID0gdGhpc1sneCddO1xuICAgIGNvbnN0IHkgPSB0aGlzWyd5J107XG4gICAgY29uc3QgeiA9IHRoaXNbJ3onXTtcblxuICAgIGxldCBub3JtID0gTWF0aC5zcXJ0KHcgKiB3ICsgeCAqIHggKyB5ICogeSArIHogKiB6KTtcblxuICAgIGlmIChub3JtIDwgRVBTSUxPTikge1xuICAgICAgcmV0dXJuIFF1YXRlcm5pb25bJ1pFUk8nXTtcbiAgICB9XG5cbiAgICBub3JtID0gMSAvIG5vcm07XG5cbiAgICByZXR1cm4gbmV3UXVhdGVybmlvbih3ICogbm9ybSwgeCAqIG5vcm0sIHkgKiBub3JtLCB6ICogbm9ybSk7XG4gIH0sXG4gIC8qKlxuICAgKiBDYWxjdWxhdGVzIHRoZSBIYW1pbHRvbiBwcm9kdWN0IG9mIHR3byBxdWF0ZXJuaW9uc1xuICAgKiBMZWF2aW5nIG91dCB0aGUgaW1hZ2luYXJ5IHBhcnQgcmVzdWx0cyBpbiBqdXN0IHNjYWxpbmcgdGhlIHF1YXRcbiAgICpcbiAgICogQHBhcmFtIHtudW1iZXJ8T2JqZWN0fHN0cmluZ30gdyByZWFsXG4gICAqIEBwYXJhbSB7bnVtYmVyPX0geCBpbWFnXG4gICAqIEBwYXJhbSB7bnVtYmVyPX0geSBpbWFnXG4gICAqIEBwYXJhbSB7bnVtYmVyPX0geiBpbWFnXG4gICAqIEByZXR1cm5zIHtRdWF0ZXJuaW9ufVxuICAgKi9cbiAgJ211bCc6IGZ1bmN0aW9uICh3LCB4LCB5LCB6KSB7XG5cbiAgICBwYXJzZShQLCB3LCB4LCB5LCB6KTtcblxuICAgIC8vIFExICogUTIgPSBbdzEgKiB3MiAtIGRvdCh2MSwgdjIpLCB3MSAqIHYyICsgdzIgKiB2MSArIGNyb3NzKHYxLCB2MildXG5cbiAgICAvLyBOb3QgY29tbXV0YXRpdmUgYmVjYXVzZSBjcm9zcyh2MSwgdjIpICE9IGNyb3NzKHYyLCB2MSkhXG5cbiAgICBjb25zdCB3MSA9IHRoaXNbJ3cnXTtcbiAgICBjb25zdCB4MSA9IHRoaXNbJ3gnXTtcbiAgICBjb25zdCB5MSA9IHRoaXNbJ3knXTtcbiAgICBjb25zdCB6MSA9IHRoaXNbJ3onXTtcblxuICAgIGNvbnN0IHcyID0gUFsndyddO1xuICAgIGNvbnN0IHgyID0gUFsneCddO1xuICAgIGNvbnN0IHkyID0gUFsneSddO1xuICAgIGNvbnN0IHoyID0gUFsneiddO1xuXG4gICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oXG4gICAgICB3MSAqIHcyIC0geDEgKiB4MiAtIHkxICogeTIgLSB6MSAqIHoyLFxuICAgICAgdzEgKiB4MiArIHgxICogdzIgKyB5MSAqIHoyIC0gejEgKiB5MixcbiAgICAgIHcxICogeTIgKyB5MSAqIHcyICsgejEgKiB4MiAtIHgxICogejIsXG4gICAgICB3MSAqIHoyICsgejEgKiB3MiArIHgxICogeTIgLSB5MSAqIHgyKTtcbiAgfSxcbiAgLyoqXG4gICAqIFNjYWxlcyBhIHF1YXRlcm5pb24gYnkgYSBzY2FsYXIsIGZhc3RlciB0aGFuIHVzaW5nIG11bHRpcGxpY2F0aW9uXG4gICAqXG4gICAqIEBwYXJhbSB7bnVtYmVyfSBzIHNjYWxpbmcgZmFjdG9yXG4gICAqIEByZXR1cm5zIHtRdWF0ZXJuaW9ufVxuICAgKi9cbiAgJ3NjYWxlJzogZnVuY3Rpb24gKHMpIHtcblxuICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKFxuICAgICAgdGhpc1sndyddICogcyxcbiAgICAgIHRoaXNbJ3gnXSAqIHMsXG4gICAgICB0aGlzWyd5J10gKiBzLFxuICAgICAgdGhpc1sneiddICogcyk7XG4gIH0sXG4gIC8qKlxuICAgKiBDYWxjdWxhdGVzIHRoZSBkb3QgcHJvZHVjdCBvZiB0d28gcXVhdGVybmlvbnNcbiAgICpcbiAgICogQHBhcmFtIHtudW1iZXJ8T2JqZWN0fHN0cmluZ30gdyByZWFsXG4gICAqIEBwYXJhbSB7bnVtYmVyPX0geCBpbWFnXG4gICAqIEBwYXJhbSB7bnVtYmVyPX0geSBpbWFnXG4gICAqIEBwYXJhbSB7bnVtYmVyPX0geiBpbWFnXG4gICAqIEByZXR1cm5zIHtudW1iZXJ9XG4gICAqL1xuICAnZG90JzogZnVuY3Rpb24gKHcsIHgsIHksIHopIHtcblxuICAgIHBhcnNlKFAsIHcsIHgsIHksIHopO1xuXG4gICAgLy8gZG90KFExLCBRMikgOj0gdzEgKiB3MiArIGRvdCh2MSwgdjIpXG5cbiAgICByZXR1cm4gdGhpc1sndyddICogUFsndyddICsgdGhpc1sneCddICogUFsneCddICsgdGhpc1sneSddICogUFsneSddICsgdGhpc1sneiddICogUFsneiddO1xuICB9LFxuICAvKipcbiAgICogQ2FsY3VsYXRlcyB0aGUgaW52ZXJzZSBvZiBhIHF1YXQgZm9yIG5vbi1ub3JtYWxpemVkIHF1YXRzIHN1Y2ggdGhhdFxuICAgKiBRXi0xICogUSA9IDEgYW5kIFEgKiBRXi0xID0gMVxuICAgKlxuICAgKiBAcmV0dXJucyB7UXVhdGVybmlvbn1cbiAgICovXG4gICdpbnZlcnNlJzogZnVuY3Rpb24gKCkge1xuXG4gICAgLy8gUV4tMSA6PSBRJyAvIHxRfF4yXG4gICAgLy8gICAgICAgPSBbdyAvICh3XjIgKyB8dnxeMiksIC12IC8gKHdeMiArIHx2fF4yKV1cblxuICAgIC8vIFByb29mOlxuICAgIC8vIFEgKiBRXi0xID0gW3csIHZdICogW3cgLyAod14yICsgfHZ8XjIpLCAtdiAvICh3XjIgKyB8dnxeMildXG4gICAgLy8gICAgICAgICAgPSBbMSwgMF1cbiAgICAvLyBRXi0xICogUSA9IFt3IC8gKHdeMiArIHx2fF4yKSwgLXYgLyAod14yICsgfHZ8XjIpXSAqIFt3LCB2XVxuICAgIC8vICAgICAgICAgID0gWzEsIDBdLlxuXG4gICAgY29uc3QgdyA9IHRoaXNbJ3cnXTtcbiAgICBjb25zdCB4ID0gdGhpc1sneCddO1xuICAgIGNvbnN0IHkgPSB0aGlzWyd5J107XG4gICAgY29uc3QgeiA9IHRoaXNbJ3onXTtcblxuICAgIGxldCBub3JtU3EgPSB3ICogdyArIHggKiB4ICsgeSAqIHkgKyB6ICogejtcblxuICAgIGlmIChub3JtU3EgPT09IDApIHtcbiAgICAgIHJldHVybiBRdWF0ZXJuaW9uWydaRVJPJ107IC8vIFRPRE86IElzIHRoZSByZXN1bHQgemVybyBvciBvbmUgd2hlbiB0aGUgbm9ybSBpcyB6ZXJvP1xuICAgIH1cblxuICAgIG5vcm1TcSA9IDEgLyBub3JtU3E7XG5cbiAgICByZXR1cm4gbmV3UXVhdGVybmlvbih3ICogbm9ybVNxLCAteCAqIG5vcm1TcSwgLXkgKiBub3JtU3EsIC16ICogbm9ybVNxKTtcbiAgfSxcbiAgLyoqXG4gICAqIE11bHRpcGxpZXMgYSBxdWF0ZXJuaW9uIHdpdGggdGhlIGludmVyc2Ugb2YgYSBzZWNvbmQgcXVhdGVybmlvblxuICAgKlxuICAgKiBAcGFyYW0ge251bWJlcnxPYmplY3R8c3RyaW5nfSB3IHJlYWxcbiAgICogQHBhcmFtIHtudW1iZXI9fSB4IGltYWdcbiAgICogQHBhcmFtIHtudW1iZXI9fSB5IGltYWdcbiAgICogQHBhcmFtIHtudW1iZXI9fSB6IGltYWdcbiAgICogQHJldHVybnMge1F1YXRlcm5pb259XG4gICAqL1xuICAnZGl2JzogZnVuY3Rpb24gKHcsIHgsIHksIHopIHtcblxuICAgIHBhcnNlKFAsIHcsIHgsIHksIHopO1xuXG4gICAgLy8gUTEgLyBRMiA6PSBRMSAqIFEyXi0xXG5cbiAgICBjb25zdCB3MSA9IHRoaXNbJ3cnXTtcbiAgICBjb25zdCB4MSA9IHRoaXNbJ3gnXTtcbiAgICBjb25zdCB5MSA9IHRoaXNbJ3knXTtcbiAgICBjb25zdCB6MSA9IHRoaXNbJ3onXTtcblxuICAgIGNvbnN0IHcyID0gUFsndyddO1xuICAgIGNvbnN0IHgyID0gUFsneCddO1xuICAgIGNvbnN0IHkyID0gUFsneSddO1xuICAgIGNvbnN0IHoyID0gUFsneiddO1xuXG4gICAgbGV0IG5vcm1TcSA9IHcyICogdzIgKyB4MiAqIHgyICsgeTIgKiB5MiArIHoyICogejI7XG5cbiAgICBpZiAobm9ybVNxID09PSAwKSB7XG4gICAgICByZXR1cm4gUXVhdGVybmlvblsnWkVSTyddOyAvLyBUT0RPOiBJcyB0aGUgcmVzdWx0IHplcm8gb3Igb25lIHdoZW4gdGhlIG5vcm0gaXMgemVybz9cbiAgICB9XG5cbiAgICBub3JtU3EgPSAxIC8gbm9ybVNxO1xuXG4gICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oXG4gICAgICAodzEgKiB3MiArIHgxICogeDIgKyB5MSAqIHkyICsgejEgKiB6MikgKiBub3JtU3EsXG4gICAgICAoeDEgKiB3MiAtIHcxICogeDIgLSB5MSAqIHoyICsgejEgKiB5MikgKiBub3JtU3EsXG4gICAgICAoeTEgKiB3MiAtIHcxICogeTIgLSB6MSAqIHgyICsgeDEgKiB6MikgKiBub3JtU3EsXG4gICAgICAoejEgKiB3MiAtIHcxICogejIgLSB4MSAqIHkyICsgeTEgKiB4MikgKiBub3JtU3EpO1xuICB9LFxuICAvKipcbiAgICogQ2FsY3VsYXRlcyB0aGUgY29uanVnYXRlIG9mIGEgcXVhdGVybmlvblxuICAgKlxuICAgKiBAcmV0dXJucyB7UXVhdGVybmlvbn1cbiAgICovXG4gICdjb25qdWdhdGUnOiBmdW5jdGlvbiAoKSB7XG5cbiAgICAvLyBRJyA9IFtzLCAtdl1cblxuICAgIC8vIElmIHRoZSBxdWF0ZXJuaW9uIGlzIG5vcm1hbGl6ZWQsXG4gICAgLy8gdGhlIGNvbmp1Z2F0ZSBpcyB0aGUgaW52ZXJzZSBvZiB0aGUgcXVhdGVybmlvbiAtIGJ1dCBmYXN0ZXJcbiAgICAvLyBRJyAqIFEgPSBRICogUScgPSAxXG5cbiAgICAvLyBBZGRpdGlvbmFsbHksIHRoZSBjb25qdWdhdGUgb2YgYSB1bml0IHF1YXRlcm5pb24gaXMgYSByb3RhdGlvbiB3aXRoIHRoZSBzYW1lXG4gICAgLy8gYW5nbGUgYnV0IHRoZSBvcHBvc2l0ZSBheGlzLlxuXG4gICAgLy8gTW9yZW92ZXIgdGhlIGZvbGxvd2luZyBwcm9wZXJ0eSBob2xkczpcbiAgICAvLyAoUTEgKiBRMiknID0gUTInICogUTEnXG5cbiAgICByZXR1cm4gbmV3UXVhdGVybmlvbih0aGlzWyd3J10sIC10aGlzWyd4J10sIC10aGlzWyd5J10sIC10aGlzWyd6J10pO1xuICB9LFxuICAvKipcbiAgICogQ2FsY3VsYXRlcyB0aGUgbmF0dXJhbCBleHBvbmVudGlhdGlvbiBvZiB0aGUgcXVhdGVybmlvblxuICAgKlxuICAgKiBAcmV0dXJucyB7UXVhdGVybmlvbn1cbiAgICovXG4gICdleHAnOiBmdW5jdGlvbiAoKSB7XG5cbiAgICBjb25zdCB3ID0gdGhpc1sndyddO1xuICAgIGNvbnN0IHggPSB0aGlzWyd4J107XG4gICAgY29uc3QgeSA9IHRoaXNbJ3knXTtcbiAgICBjb25zdCB6ID0gdGhpc1sneiddO1xuXG4gICAgY29uc3Qgdk5vcm0gPSBNYXRoLnNxcnQoeCAqIHggKyB5ICogeSArIHogKiB6KTtcbiAgICBjb25zdCB3RXhwID0gTWF0aC5leHAodyk7XG4gICAgY29uc3Qgc2NhbGUgPSB3RXhwICogTWF0aC5zaW4odk5vcm0pIC8gdk5vcm07XG5cbiAgICBpZiAodk5vcm0gPT09IDApIHtcbiAgICAgIC8vcmV0dXJuIG5ld1F1YXRlcm5pb24od0V4cCAqIE1hdGguY29zKHZOb3JtKSwgMCwgMCwgMCk7XG4gICAgICByZXR1cm4gbmV3UXVhdGVybmlvbih3RXhwLCAwLCAwLCAwKTtcbiAgICB9XG5cbiAgICByZXR1cm4gbmV3UXVhdGVybmlvbihcbiAgICAgIHdFeHAgKiBNYXRoLmNvcyh2Tm9ybSksXG4gICAgICB4ICogc2NhbGUsXG4gICAgICB5ICogc2NhbGUsXG4gICAgICB6ICogc2NhbGUpO1xuICB9LFxuICAvKipcbiAgICogQ2FsY3VsYXRlcyB0aGUgbmF0dXJhbCBsb2dhcml0aG0gb2YgdGhlIHF1YXRlcm5pb25cbiAgICpcbiAgICogQHJldHVybnMge1F1YXRlcm5pb259XG4gICAqL1xuICAnbG9nJzogZnVuY3Rpb24gKCkge1xuXG4gICAgY29uc3QgdyA9IHRoaXNbJ3cnXTtcbiAgICBjb25zdCB4ID0gdGhpc1sneCddO1xuICAgIGNvbnN0IHkgPSB0aGlzWyd5J107XG4gICAgY29uc3QgeiA9IHRoaXNbJ3onXTtcblxuICAgIGlmICh5ID09PSAwICYmIHogPT09IDApIHtcbiAgICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKFxuICAgICAgICBsb2dIeXBvdCh3LCB4KSxcbiAgICAgICAgTWF0aC5hdGFuMih4LCB3KSwgMCwgMCk7XG4gICAgfVxuXG4gICAgY29uc3QgcU5vcm0yID0geCAqIHggKyB5ICogeSArIHogKiB6ICsgdyAqIHc7XG4gICAgY29uc3Qgdk5vcm0gPSBNYXRoLnNxcnQoeCAqIHggKyB5ICogeSArIHogKiB6KTtcblxuICAgIGNvbnN0IHNjYWxlID0gTWF0aC5hdGFuMih2Tm9ybSwgdykgLyB2Tm9ybTsgLy8gQWx0ZXJuYXRpdmU6IGFjb3ModyAvIHFOb3JtKSAvIHZOb3JtXG5cbiAgICByZXR1cm4gbmV3UXVhdGVybmlvbihcbiAgICAgIE1hdGgubG9nKHFOb3JtMikgKiAwLjUsXG4gICAgICB4ICogc2NhbGUsXG4gICAgICB5ICogc2NhbGUsXG4gICAgICB6ICogc2NhbGUpO1xuICB9LFxuICAvKipcbiAgICogQ2FsY3VsYXRlcyB0aGUgcG93ZXIgb2YgYSBxdWF0ZXJuaW9uIHJhaXNlZCB0byBhIHJlYWwgbnVtYmVyIG9yIGFub3RoZXIgcXVhdGVybmlvblxuICAgKlxuICAgKiBAcGFyYW0ge251bWJlcnxPYmplY3R8c3RyaW5nfSB3IHJlYWxcbiAgICogQHBhcmFtIHtudW1iZXI9fSB4IGltYWdcbiAgICogQHBhcmFtIHtudW1iZXI9fSB5IGltYWdcbiAgICogQHBhcmFtIHtudW1iZXI9fSB6IGltYWdcbiAgICogQHJldHVybnMge1F1YXRlcm5pb259XG4gICAqL1xuICAncG93JzogZnVuY3Rpb24gKHcsIHgsIHksIHopIHtcblxuICAgIHBhcnNlKFAsIHcsIHgsIHksIHopO1xuXG4gICAgaWYgKFBbJ3knXSA9PT0gMCAmJiBQWyd6J10gPT09IDApIHtcblxuICAgICAgaWYgKFBbJ3cnXSA9PT0gMSAmJiBQWyd4J10gPT09IDApIHtcbiAgICAgICAgcmV0dXJuIHRoaXM7XG4gICAgICB9XG5cbiAgICAgIGlmIChQWyd3J10gPT09IDAgJiYgUFsneCddID09PSAwKSB7XG4gICAgICAgIHJldHVybiBRdWF0ZXJuaW9uWydPTkUnXTtcbiAgICAgIH1cblxuICAgICAgLy8gQ2hlY2sgaWYgd2UgY2FuIG9wZXJhdGUgaW4gQ1xuICAgICAgLy8gQm9ycm93ZWQgZnJvbSBjb21wbGV4LmpzXG4gICAgICBpZiAodGhpc1sneSddID09PSAwICYmIHRoaXNbJ3onXSA9PT0gMCkge1xuXG4gICAgICAgIGxldCBhID0gdGhpc1sndyddO1xuICAgICAgICBsZXQgYiA9IHRoaXNbJ3gnXTtcblxuICAgICAgICBpZiAoYSA9PT0gMCAmJiBiID09PSAwKSB7XG4gICAgICAgICAgcmV0dXJuIFF1YXRlcm5pb25bJ1pFUk8nXTtcbiAgICAgICAgfVxuXG4gICAgICAgIGxldCBhcmcgPSBNYXRoLmF0YW4yKGIsIGEpO1xuICAgICAgICBsZXQgbG9oID0gbG9nSHlwb3QoYSwgYik7XG5cbiAgICAgICAgaWYgKFBbJ3gnXSA9PT0gMCkge1xuXG4gICAgICAgICAgaWYgKGIgPT09IDAgJiYgYSA+PSAwKSB7XG5cbiAgICAgICAgICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKE1hdGgucG93KGEsIFBbJ3cnXSksIDAsIDAsIDApO1xuXG4gICAgICAgICAgfSBlbHNlIGlmIChhID09PSAwKSB7XG5cbiAgICAgICAgICAgIHN3aXRjaCAoUFsndyddICUgNCkge1xuICAgICAgICAgICAgICBjYXNlIDA6XG4gICAgICAgICAgICAgICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oTWF0aC5wb3coYiwgUFsndyddKSwgMCwgMCwgMCk7XG4gICAgICAgICAgICAgIGNhc2UgMTpcbiAgICAgICAgICAgICAgICByZXR1cm4gbmV3UXVhdGVybmlvbigwLCBNYXRoLnBvdyhiLCBQWyd3J10pLCAwLCAwKTtcbiAgICAgICAgICAgICAgY2FzZSAyOlxuICAgICAgICAgICAgICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKC1NYXRoLnBvdyhiLCBQWyd3J10pLCAwLCAwLCAwKTtcbiAgICAgICAgICAgICAgY2FzZSAzOlxuICAgICAgICAgICAgICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKDAsIC1NYXRoLnBvdyhiLCBQWyd3J10pLCAwLCAwKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICBhID0gTWF0aC5leHAoUFsndyddICogbG9oIC0gUFsneCddICogYXJnKTtcbiAgICAgICAgYiA9IFBbJ3gnXSAqIGxvaCArIFBbJ3cnXSAqIGFyZztcbiAgICAgICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oXG4gICAgICAgICAgYSAqIE1hdGguY29zKGIpLFxuICAgICAgICAgIGEgKiBNYXRoLnNpbihiKSwgMCwgMCk7XG4gICAgICB9XG4gICAgfVxuXG4gICAgLy8gTm9ybWFsIHF1YXRlcm5pb24gYmVoYXZpb3JcbiAgICAvLyBxXnAgPSBlXmxuKHFecCkgPSBlXihsbihxKSpwKVxuICAgIHJldHVybiB0aGlzWydsb2cnXSgpWydtdWwnXShQKVsnZXhwJ10oKTtcbiAgfSxcbiAgLyoqXG4gICAqIENoZWNrcyBpZiB0d28gcXVhdHMgYXJlIHRoZSBzYW1lXG4gICAqXG4gICAqIEBwYXJhbSB7bnVtYmVyfE9iamVjdHxzdHJpbmd9IHcgcmVhbFxuICAgKiBAcGFyYW0ge251bWJlcj19IHggaW1hZ1xuICAgKiBAcGFyYW0ge251bWJlcj19IHkgaW1hZ1xuICAgKiBAcGFyYW0ge251bWJlcj19IHogaW1hZ1xuICAgKiBAcmV0dXJucyB7Ym9vbGVhbn1cbiAgICovXG4gICdlcXVhbHMnOiBmdW5jdGlvbiAodywgeCwgeSwgeikge1xuXG4gICAgcGFyc2UoUCwgdywgeCwgeSwgeik7XG5cbiAgICBjb25zdCBlcHMgPSBFUFNJTE9OO1xuXG4gICAgLy8gbWF5YmUgY2hlY2sgZm9yIE5hTidzIGhlcmU/XG4gICAgcmV0dXJuIE1hdGguYWJzKFBbJ3cnXSAtIHRoaXNbJ3cnXSkgPCBlcHNcbiAgICAgICYmIE1hdGguYWJzKFBbJ3gnXSAtIHRoaXNbJ3gnXSkgPCBlcHNcbiAgICAgICYmIE1hdGguYWJzKFBbJ3knXSAtIHRoaXNbJ3knXSkgPCBlcHNcbiAgICAgICYmIE1hdGguYWJzKFBbJ3onXSAtIHRoaXNbJ3onXSkgPCBlcHM7XG4gIH0sXG4gIC8qKlxuICAgKiBDaGVja3MgaWYgYWxsIHBhcnRzIG9mIGEgcXVhdGVybmlvbiBhcmUgZmluaXRlXG4gICAqXG4gICAqIEByZXR1cm5zIHtib29sZWFufVxuICAgKi9cbiAgJ2lzRmluaXRlJzogZnVuY3Rpb24gKCkge1xuXG4gICAgcmV0dXJuIGlzRmluaXRlKHRoaXNbJ3cnXSkgJiYgaXNGaW5pdGUodGhpc1sneCddKSAmJiBpc0Zpbml0ZSh0aGlzWyd5J10pICYmIGlzRmluaXRlKHRoaXNbJ3onXSk7XG4gIH0sXG4gIC8qKlxuICAgKiBDaGVja3MgaWYgYW55IG9mIHRoZSBwYXJ0cyBvZiB0aGUgcXVhdGVybmlvbiBpcyBub3QgYSBudW1iZXJcbiAgICpcbiAgICogQHJldHVybnMge2Jvb2xlYW59XG4gICAqL1xuICAnaXNOYU4nOiBmdW5jdGlvbiAoKSB7XG5cbiAgICByZXR1cm4gaXNOYU4odGhpc1sndyddKSB8fCBpc05hTih0aGlzWyd4J10pIHx8IGlzTmFOKHRoaXNbJ3knXSkgfHwgaXNOYU4odGhpc1sneiddKTtcbiAgfSxcbiAgLyoqXG4gICAqIEdldHMgdGhlIFF1YXRlcm5pb24gYXMgYSB3ZWxsIGZvcm1hdHRlZCBzdHJpbmdcbiAgICpcbiAgICogQHJldHVybnMge3N0cmluZ31cbiAgICovXG4gICd0b1N0cmluZyc6IGZ1bmN0aW9uICgpIHtcblxuICAgIGNvbnN0IHcgPSB0aGlzWyd3J107XG4gICAgY29uc3QgeCA9IHRoaXNbJ3gnXTtcbiAgICBjb25zdCB5ID0gdGhpc1sneSddO1xuICAgIGNvbnN0IHogPSB0aGlzWyd6J107XG4gICAgbGV0IHJldCA9ICcnO1xuXG4gICAgaWYgKGlzTmFOKHcpIHx8IGlzTmFOKHgpIHx8IGlzTmFOKHkpIHx8IGlzTmFOKHopKSB7XG4gICAgICByZXR1cm4gJ05hTic7XG4gICAgfVxuXG4gICAgLy8gQWx0ZXJuYXRpdmUgZGVzaWduP1xuICAgIC8vICcoJWYsIFslZiAlZiAlZl0pJ1xuXG4gICAgcmV0ID0gbnVtVG9TdHIodywgJycsIHJldCk7XG4gICAgcmV0ICs9IG51bVRvU3RyKHgsICdpJywgcmV0KTtcbiAgICByZXQgKz0gbnVtVG9TdHIoeSwgJ2onLCByZXQpO1xuICAgIHJldCArPSBudW1Ub1N0cih6LCAnaycsIHJldCk7XG5cbiAgICBpZiAoJycgPT09IHJldClcbiAgICAgIHJldHVybiAnMCc7XG5cbiAgICByZXR1cm4gcmV0O1xuICB9LFxuICAvKipcbiAgICogUmV0dXJucyB0aGUgcmVhbCBwYXJ0IG9mIHRoZSBxdWF0ZXJuaW9uXG4gICAqXG4gICAqIEByZXR1cm5zIHtudW1iZXJ9XG4gICAqL1xuICAncmVhbCc6IGZ1bmN0aW9uICgpIHtcblxuICAgIHJldHVybiB0aGlzWyd3J107XG4gIH0sXG4gIC8qKlxuICAgKiBSZXR1cm5zIHRoZSBpbWFnaW5hcnkgcGFydCBvZiB0aGUgcXVhdGVybmlvbiBhcyBhIDNEIHZlY3RvciAvIGFycmF5XG4gICAqXG4gICAqIEByZXR1cm5zIHtBcnJheX1cbiAgICovXG4gICdpbWFnJzogZnVuY3Rpb24gKCkge1xuXG4gICAgcmV0dXJuIFt0aGlzWyd4J10sIHRoaXNbJ3knXSwgdGhpc1sneiddXTtcbiAgfSxcbiAgLyoqXG4gICAqIEdldHMgdGhlIGFjdHVhbCBxdWF0ZXJuaW9uIGFzIGEgNEQgdmVjdG9yIC8gYXJyYXlcbiAgICpcbiAgICogQHJldHVybnMge0FycmF5fVxuICAgKi9cbiAgJ3RvVmVjdG9yJzogZnVuY3Rpb24gKCkge1xuXG4gICAgcmV0dXJuIFt0aGlzWyd3J10sIHRoaXNbJ3gnXSwgdGhpc1sneSddLCB0aGlzWyd6J11dO1xuICB9LFxuICAvKipcbiAgICogQ2FsY3VsYXRlcyB0aGUgM3gzIHJvdGF0aW9uIG1hdHJpeCBmb3IgdGhlIGN1cnJlbnQgcXVhdFxuICAgKlxuICAgKiBAcGFyYW0ge2Jvb2xlYW49fSB0d29EXG4gICAqIEByZXR1cm5zIHtBcnJheX1cbiAgICovXG4gICd0b01hdHJpeCc6IGZ1bmN0aW9uICh0d29EKSB7XG5cbiAgICBjb25zdCB3ID0gdGhpc1sndyddO1xuICAgIGNvbnN0IHggPSB0aGlzWyd4J107XG4gICAgY29uc3QgeSA9IHRoaXNbJ3knXTtcbiAgICBjb25zdCB6ID0gdGhpc1sneiddO1xuXG4gICAgY29uc3Qgd3ggPSB3ICogeCwgd3kgPSB3ICogeSwgd3ogPSB3ICogejtcbiAgICBjb25zdCB4eCA9IHggKiB4LCB4eSA9IHggKiB5LCB4eiA9IHggKiB6O1xuICAgIGNvbnN0IHl5ID0geSAqIHksIHl6ID0geSAqIHosIHp6ID0geiAqIHo7XG5cbiAgICBpZiAodHdvRCkge1xuICAgICAgcmV0dXJuIFtcbiAgICAgICAgWzEgLSAyICogKHl5ICsgenopLCAyICogKHh5IC0gd3opLCAyICogKHh6ICsgd3kpXSxcbiAgICAgICAgWzIgKiAoeHkgKyB3eiksIDEgLSAyICogKHh4ICsgenopLCAyICogKHl6IC0gd3gpXSxcbiAgICAgICAgWzIgKiAoeHogLSB3eSksIDIgKiAoeXogKyB3eCksIDEgLSAyICogKHh4ICsgeXkpXV07XG4gICAgfVxuXG4gICAgcmV0dXJuIFtcbiAgICAgIDEgLSAyICogKHl5ICsgenopLCAyICogKHh5IC0gd3opLCAyICogKHh6ICsgd3kpLFxuICAgICAgMiAqICh4eSArIHd6KSwgMSAtIDIgKiAoeHggKyB6eiksIDIgKiAoeXogLSB3eCksXG4gICAgICAyICogKHh6IC0gd3kpLCAyICogKHl6ICsgd3gpLCAxIC0gMiAqICh4eCArIHl5KV07XG4gIH0sXG4gIC8qKlxuICAgKiBDYWxjdWxhdGVzIHRoZSBob21vZ2VuZW91cyA0eDQgcm90YXRpb24gbWF0cml4IGZvciB0aGUgY3VycmVudCBxdWF0XG4gICAqXG4gICAqIEBwYXJhbSB7Ym9vbGVhbj19IHR3b0RcbiAgICogQHJldHVybnMge0FycmF5fVxuICAgKi9cbiAgJ3RvTWF0cml4NCc6IGZ1bmN0aW9uICh0d29EKSB7XG5cbiAgICBjb25zdCB3ID0gdGhpc1sndyddO1xuICAgIGNvbnN0IHggPSB0aGlzWyd4J107XG4gICAgY29uc3QgeSA9IHRoaXNbJ3knXTtcbiAgICBjb25zdCB6ID0gdGhpc1sneiddO1xuXG4gICAgY29uc3Qgd3ggPSB3ICogeCwgd3kgPSB3ICogeSwgd3ogPSB3ICogejtcbiAgICBjb25zdCB4eCA9IHggKiB4LCB4eSA9IHggKiB5LCB4eiA9IHggKiB6O1xuICAgIGNvbnN0IHl5ID0geSAqIHksIHl6ID0geSAqIHosIHp6ID0geiAqIHo7XG5cbiAgICBpZiAodHdvRCkge1xuICAgICAgcmV0dXJuIFtcbiAgICAgICAgWzEgLSAyICogKHl5ICsgenopLCAyICogKHh5IC0gd3opLCAyICogKHh6ICsgd3kpLCAwXSxcbiAgICAgICAgWzIgKiAoeHkgKyB3eiksIDEgLSAyICogKHh4ICsgenopLCAyICogKHl6IC0gd3gpLCAwXSxcbiAgICAgICAgWzIgKiAoeHogLSB3eSksIDIgKiAoeXogKyB3eCksIDEgLSAyICogKHh4ICsgeXkpLCAwXSxcbiAgICAgICAgWzAsIDAsIDAsIDFdXTtcbiAgICB9XG5cbiAgICByZXR1cm4gW1xuICAgICAgMSAtIDIgKiAoeXkgKyB6eiksIDIgKiAoeHkgLSB3eiksIDIgKiAoeHogKyB3eSksIDAsXG4gICAgICAyICogKHh5ICsgd3opLCAxIC0gMiAqICh4eCArIHp6KSwgMiAqICh5eiAtIHd4KSwgMCxcbiAgICAgIDIgKiAoeHogLSB3eSksIDIgKiAoeXogKyB3eCksIDEgLSAyICogKHh4ICsgeXkpLCAwLFxuICAgICAgMCwgMCwgMCwgMV07XG4gIH0sXG4gIC8qKlxuICAgKiBEZXRlcm1pbmVzIHRoZSBob21vZ2VuZW91cyByb3RhdGlvbiBtYXRyaXggc3RyaW5nIHVzZWQgZm9yIENTUyAzRCB0cmFuc2Zvcm1zXG4gICAqXG4gICAqIEByZXR1cm5zIHtzdHJpbmd9XG4gICAqL1xuICAndG9DU1NUcmFuc2Zvcm0nOiBmdW5jdGlvbiAoKSB7XG5cbiAgICBjb25zdCB3ID0gdGhpc1sndyddO1xuXG4gICAgbGV0IGFuZ2xlID0gMiAqIE1hdGguYWNvcyh3KTtcbiAgICBsZXQgc2luMiA9IDEgLSB3ICogdztcblxuICAgIGlmIChzaW4yIDwgRVBTSUxPTikge1xuICAgICAgYW5nbGUgPSAwO1xuICAgICAgc2luMiA9IDE7XG4gICAgfSBlbHNlIHtcbiAgICAgIHNpbjIgPSAxIC8gTWF0aC5zcXJ0KHNpbjIpOyAvLyBSZS11c2UgdmFyaWFibGUgc2luXjIgZm9yIDEgLyBzaW5cbiAgICB9XG4gICAgcmV0dXJuIFwicm90YXRlM2QoXCIgKyB0aGlzWyd4J10gKiBzaW4yICsgXCIsXCIgKyB0aGlzWyd5J10gKiBzaW4yICsgXCIsXCIgKyB0aGlzWyd6J10gKiBzaW4yICsgXCIsXCIgKyBhbmdsZSArIFwicmFkKVwiO1xuICB9LFxuICAvKipcbiAgICogQ2FsY3VsYXRlcyB0aGUgYXhpcyBhbmQgYW5nbGUgcmVwcmVzZW50YXRpb24gb2YgdGhlIHF1YXRlcm5pb25cbiAgICpcbiAgICogQHJldHVybnMge0FycmF5fVxuICAgKi9cbiAgJ3RvQXhpc0FuZ2xlJzogZnVuY3Rpb24gKCkge1xuXG4gICAgY29uc3QgdyA9IHRoaXNbJ3cnXTtcbiAgICBjb25zdCBzaW4yID0gMSAtIHcgKiB3OyAvLyBzaW4oYW5nbGUgLyAyKSA9IHNpbihhY29zKHcpKSA9IHNxcnQoMSAtIHdeMikgPSB8dnwsIHNpbmNlIDEgPSBkb3QoUSkgPD0+IGRvdCh2KSA9IDEgLSB3XjJcblxuICAgIGlmIChzaW4yIDwgRVBTSUxPTikgeyAvLyBBbHRlcm5hdGl2ZWx5IHx2fCA9PSAwXG4gICAgICAvLyBJZiB0aGUgc2luZSBpcyBjbG9zZSB0byAwLCB3ZSdyZSBjbG9zZSB0byB0aGUgdW5pdCBxdWF0ZXJuaW9uIGFuZCB0aGUgZGlyZWN0aW9uIGRvZXMgbm90IG1hdHRlclxuICAgICAgcmV0dXJuIFtbdGhpc1sneCddLCB0aGlzWyd5J10sIHRoaXNbJ3onXV0sIDBdOyAvLyBvciBbWzEsIDAsIDBdLCAwXSA/ICBvciBbWzAsIDAsIDBdLCAwXSA/XG4gICAgfVxuXG4gICAgY29uc3QgaXNpbiA9IDEgLyBNYXRoLnNxcnQoc2luMik7XG4gICAgY29uc3QgYW5nbGUgPSAyICogTWF0aC5hY29zKHcpOyAvLyBBbHRlcm5hdGl2ZWx5OiAyICogYXRhbjIofHZ8LCB3KVxuICAgIHJldHVybiBbW3RoaXNbJ3gnXSAqIGlzaW4sIHRoaXNbJ3knXSAqIGlzaW4sIHRoaXNbJ3onXSAqIGlzaW5dLCBhbmdsZV07XG4gIH0sXG4gIC8qKlxuICAgKiBDYWxjdWxhdGVzIHRoZSBFdWxlciBhbmdsZXMgcmVwcmVzZW50ZWQgYnkgdGhlIGN1cnJlbnQgcXVhdCAobXVsdGlwbGljYXRpb24gb3JkZXIgZnJvbSByaWdodCB0byBsZWZ0KVxuICAgKiBcbiAgICogQHBhcmFtIHtzdHJpbmc9fSBvcmRlciBBeGlzIG9yZGVyIChUYWl0IEJyeWFuKVxuICAgKiBAcmV0dXJucyB7T2JqZWN0fVxuICAgKi9cbiAgJ3RvRXVsZXInOiBmdW5jdGlvbiAob3JkZXIpIHtcblxuICAgIGNvbnN0IHcgPSB0aGlzWyd3J107XG4gICAgY29uc3QgeCA9IHRoaXNbJ3gnXTtcbiAgICBjb25zdCB5ID0gdGhpc1sneSddO1xuICAgIGNvbnN0IHogPSB0aGlzWyd6J107XG5cbiAgICBjb25zdCB3eCA9IHcgKiB4LCB3eSA9IHcgKiB5LCB3eiA9IHcgKiB6O1xuICAgIGNvbnN0IHh4ID0geCAqIHgsIHh5ID0geCAqIHksIHh6ID0geCAqIHo7XG4gICAgY29uc3QgeXkgPSB5ICogeSwgeXogPSB5ICogeiwgenogPSB6ICogejtcblxuICAgIGZ1bmN0aW9uIGFzaW4odCkge1xuICAgICAgcmV0dXJuIHQgPj0gMSA/IE1hdGguUEkgLyAyIDogKHQgPD0gLTEgPyAtTWF0aC5QSSAvIDIgOiBNYXRoLmFzaW4odCkpO1xuICAgIH1cblxuICAgIGlmIChvcmRlciA9PT0gdW5kZWZpbmVkIHx8IG9yZGVyID09PSAnWlhZJykge1xuICAgICAgcmV0dXJuIFtcbiAgICAgICAgLU1hdGguYXRhbjIoMiAqICh4eSAtIHd6KSwgMSAtIDIgKiAoeHggKyB6eikpLFxuICAgICAgICBhc2luKDIgKiAoeXogKyB3eCkpLFxuICAgICAgICAtTWF0aC5hdGFuMigyICogKHh6IC0gd3kpLCAxIC0gMiAqICh4eCArIHl5KSksXG4gICAgICBdO1xuICAgIH1cblxuICAgIGlmIChvcmRlciA9PT0gJ1hZWicgfHwgb3JkZXIgPT09ICdSUFknKSB7XG4gICAgICByZXR1cm4gW1xuICAgICAgICAtTWF0aC5hdGFuMigyICogKHl6IC0gd3gpLCAxIC0gMiAqICh4eCArIHl5KSksXG4gICAgICAgIGFzaW4oMiAqICh4eiArIHd5KSksXG4gICAgICAgIC1NYXRoLmF0YW4yKDIgKiAoeHkgLSB3eiksIDEgLSAyICogKHl5ICsgenopKSxcbiAgICAgIF07XG4gICAgfVxuXG4gICAgaWYgKG9yZGVyID09PSAnWVhaJykge1xuICAgICAgcmV0dXJuIFtcbiAgICAgICAgTWF0aC5hdGFuMigyICogKHh6ICsgd3kpLCAxIC0gMiAqICh4eCArIHl5KSksXG4gICAgICAgIC1hc2luKDIgKiAoeXogLSB3eCkpLFxuICAgICAgICBNYXRoLmF0YW4yKDIgKiAoeHkgKyB3eiksIDEgLSAyICogKHh4ICsgenopKSxcbiAgICAgIF07XG4gICAgfVxuXG4gICAgaWYgKG9yZGVyID09PSAnWllYJyB8fCBvcmRlciA9PT0gJ1lQUicpIHsgIC8vIHJvbGwgYXJvdW5kIFgsIHBpdGNoIGFyb3VuZCBZLCB5YXcgYXJvdW5kIFpcbiAgICAgIC8qXG4gICAgICBpZiAoMiAqICh4eiAtIHd5KSA+IC45OTkpIHtcbiAgICAgICAgcmV0dXJuIFtcbiAgICAgICAgICAyICogTWF0aC5hdGFuMih4LCB3KSxcbiAgICAgICAgICAtTWF0aC5QSSAvIDIsXG4gICAgICAgICAgMFxuICAgICAgICBdO1xuICAgICAgfVxuXG4gICAgICBpZiAoMiAqICh4eiAtIHd5KSA8IC0uOTk5KSB7XG4gICAgICAgIHJldHVybiBbXG4gICAgICAgICAgLTIgKiBNYXRoLmF0YW4yKHgsIHcpLFxuICAgICAgICAgIE1hdGguUEkgLyAyLFxuICAgICAgICAgIDBcbiAgICAgICAgXTtcbiAgICAgIH1cbiAgICAgICovXG4gICAgICByZXR1cm4gW1xuICAgICAgICBNYXRoLmF0YW4yKDIgKiAoeHkgKyB3eiksIDEgLSAyICogKHl5ICsgenopKSwgLy8gSGVhZGluZyAvIFlhd1xuICAgICAgICAtYXNpbigyICogKHh6IC0gd3kpKSwgLy8gQXR0aXR1ZGUgLyBQaXRjaFxuICAgICAgICBNYXRoLmF0YW4yKDIgKiAoeXogKyB3eCksIDEgLSAyICogKHh4ICsgeXkpKSwgLy8gQmFuayAvIFJvbGxcbiAgICAgIF07XG4gICAgfVxuXG4gICAgaWYgKG9yZGVyID09PSAnWVpYJykge1xuICAgICAgLypcbiAgICAgIGlmICgyICogKHh5ICsgd3opID4gLjk5OSkgeyAvLyBOb3J0aCBwb2xlXG4gICAgICAgIHJldHVybiBbXG4gICAgICAgICAgMiAqIE1hdGguYXRhbjIoeCwgdyksXG4gICAgICAgICAgTWF0aC5QSSAvIDIsXG4gICAgICAgICAgMFxuICAgICAgICBdO1xuICAgICAgfVxuXG4gICAgICBpZiAoMiAqICh4eSArIHd6KSA8IC0uOTk5KSB7IC8vIFNvdXRoIHBvbGVcbiAgICAgICAgcmV0dXJuIFtcbiAgICAgICAgICAtMiAqIE1hdGguYXRhbjIoeCwgdyksXG4gICAgICAgICAgLU1hdGguUEkgLyAyLFxuICAgICAgICAgIDBcbiAgICAgICAgXTtcbiAgICAgIH1cbiAgICAgICovXG4gICAgICByZXR1cm4gW1xuICAgICAgICAtTWF0aC5hdGFuMigyICogKHh6IC0gd3kpLCAxIC0gMiAqICh5eSArIHp6KSksIC8vIEhlYWRpbmdcbiAgICAgICAgYXNpbigyICogKHh5ICsgd3opKSwgLy8gQXR0aXR1ZGVcbiAgICAgICAgLU1hdGguYXRhbjIoMiAqICh5eiAtIHd4KSwgMSAtIDIgKiAoeHggKyB6eikpLCAvLyBCYW5rXG4gICAgICBdO1xuICAgIH1cblxuICAgIGlmIChvcmRlciA9PT0gJ1haWScpIHtcbiAgICAgIHJldHVybiBbXG4gICAgICAgIE1hdGguYXRhbjIoMiAqICh5eiArIHd4KSwgMSAtIDIgKiAoeHggKyB6eikpLFxuICAgICAgICAtYXNpbigyICogKHh5IC0gd3opKSxcbiAgICAgICAgTWF0aC5hdGFuMigyICogKHh6ICsgd3kpLCAxIC0gMiAqICh5eSArIHp6KSksXG4gICAgICBdO1xuICAgIH1cbiAgICByZXR1cm4gbnVsbDtcbiAgfSxcbiAgLyoqXG4gICAqIENsb25lcyB0aGUgYWN0dWFsIG9iamVjdFxuICAgKlxuICAgKiBAcmV0dXJucyB7UXVhdGVybmlvbn1cbiAgICovXG4gICdjbG9uZSc6IGZ1bmN0aW9uICgpIHtcblxuICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKHRoaXNbJ3cnXSwgdGhpc1sneCddLCB0aGlzWyd5J10sIHRoaXNbJ3onXSk7XG4gIH0sXG4gIC8qKlxuICAgKiBSb3RhdGVzIGEgdmVjdG9yIGFjY29yZGluZyB0byB0aGUgY3VycmVudCBxdWF0ZXJuaW9uLCBhc3N1bWVzIHxxfD0xXG4gICAqIEBsaW5rIGh0dHBzOi8vcmF3Lm9yZy9wcm9vZi92ZWN0b3Itcm90YXRpb24tdXNpbmctcXVhdGVybmlvbnMvXG4gICAqXG4gICAqIEBwYXJhbSB7QXJyYXl9IHYgVGhlIHZlY3RvciB0byBiZSByb3RhdGVkXG4gICAqIEByZXR1cm5zIHtBcnJheX1cbiAgICovXG4gICdyb3RhdGVWZWN0b3InOiBmdW5jdGlvbiAodikge1xuXG4gICAgY29uc3QgcXcgPSB0aGlzWyd3J107XG4gICAgY29uc3QgcXggPSB0aGlzWyd4J107XG4gICAgY29uc3QgcXkgPSB0aGlzWyd5J107XG4gICAgY29uc3QgcXogPSB0aGlzWyd6J107XG5cbiAgICBjb25zdCB2eCA9IHZbMF07XG4gICAgY29uc3QgdnkgPSB2WzFdO1xuICAgIGNvbnN0IHZ6ID0gdlsyXTtcblxuICAgIC8vIHQgPSBxIHggdlxuICAgIGxldCB0eCA9IHF5ICogdnogLSBxeiAqIHZ5O1xuICAgIGxldCB0eSA9IHF6ICogdnggLSBxeCAqIHZ6O1xuICAgIGxldCB0eiA9IHF4ICogdnkgLSBxeSAqIHZ4O1xuXG4gICAgLy8gdCA9IDJ0XG4gICAgdHggPSB0eCArIHR4O1xuICAgIHR5ID0gdHkgKyB0eTtcbiAgICB0eiA9IHR6ICsgdHo7XG5cbiAgICAvLyB2ICsgdyB0ICsgcSB4IHRcbiAgICByZXR1cm4gW1xuICAgICAgdnggKyBxdyAqIHR4ICsgcXkgKiB0eiAtIHF6ICogdHksXG4gICAgICB2eSArIHF3ICogdHkgKyBxeiAqIHR4IC0gcXggKiB0eixcbiAgICAgIHZ6ICsgcXcgKiB0eiArIHF4ICogdHkgLSBxeSAqIHR4XTtcbiAgfSxcblxuICAvKipcbiAgICogR2V0cyBhIGZ1bmN0aW9uIHRvIHNwaGVyaWNhbGx5IGludGVycG9sYXRlIGJldHdlZW4gdHdvIHF1YXRlcm5pb25zXG4gICAqIFxuICAgKiBAcmV0dXJucyBGdW5jdGlvblxuICAgKi9cbiAgJ3NsZXJwJzogZnVuY3Rpb24gKHcsIHgsIHksIHopIHtcblxuICAgIHBhcnNlKFAsIHcsIHgsIHksIHopO1xuXG4gICAgLy8gc2xlcnAoUTEsIFEyLCB0KSA6PSBRMShRMV4tMSBRMiledFxuXG4gICAgbGV0IHcxID0gdGhpc1sndyddO1xuICAgIGxldCB4MSA9IHRoaXNbJ3gnXTtcbiAgICBsZXQgeTEgPSB0aGlzWyd5J107XG4gICAgbGV0IHoxID0gdGhpc1sneiddO1xuXG4gICAgbGV0IHcyID0gUFsndyddO1xuICAgIGxldCB4MiA9IFBbJ3gnXTtcbiAgICBsZXQgeTIgPSBQWyd5J107XG4gICAgbGV0IHoyID0gUFsneiddO1xuXG4gICAgbGV0IGNvc1RoZXRhMCA9IHcxICogdzIgKyB4MSAqIHgyICsgeTEgKiB5MiArIHoxICogejI7XG5cbiAgICBpZiAoY29zVGhldGEwIDwgMCkge1xuICAgICAgdzEgPSAtdzE7XG4gICAgICB4MSA9IC14MTtcbiAgICAgIHkxID0gLXkxO1xuICAgICAgejEgPSAtejE7XG4gICAgICBjb3NUaGV0YTAgPSAtY29zVGhldGEwO1xuICAgIH1cblxuICAgIGlmIChjb3NUaGV0YTAgPj0gMSAtIEVQU0lMT04pIHtcbiAgICAgIHJldHVybiBmdW5jdGlvbiAocGN0KSB7XG4gICAgICAgIHJldHVybiBuZXdOb3JtYWxpemVkKFxuICAgICAgICAgIHcxICsgcGN0ICogKHcyIC0gdzEpLFxuICAgICAgICAgIHgxICsgcGN0ICogKHgyIC0geDEpLFxuICAgICAgICAgIHkxICsgcGN0ICogKHkyIC0geTEpLFxuICAgICAgICAgIHoxICsgcGN0ICogKHoyIC0gejEpKTtcbiAgICAgIH07XG4gICAgfVxuXG4gICAgbGV0IFRoZXRhMCA9IE1hdGguYWNvcyhjb3NUaGV0YTApO1xuICAgIGxldCBzaW5UaGV0YTAgPSBNYXRoLnNpbihUaGV0YTApO1xuXG4gICAgcmV0dXJuIGZ1bmN0aW9uIChwY3QpIHtcblxuICAgICAgbGV0IFRoZXRhID0gVGhldGEwICogcGN0O1xuICAgICAgbGV0IHNpblRoZXRhID0gTWF0aC5zaW4oVGhldGEpO1xuICAgICAgbGV0IGNvc1RoZXRhID0gTWF0aC5jb3MoVGhldGEpO1xuXG4gICAgICBsZXQgczAgPSBjb3NUaGV0YSAtIGNvc1RoZXRhMCAqIHNpblRoZXRhIC8gc2luVGhldGEwO1xuICAgICAgbGV0IHMxID0gc2luVGhldGEgLyBzaW5UaGV0YTA7XG5cbiAgICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKFxuICAgICAgICBzMCAqIHcxICsgczEgKiB3MixcbiAgICAgICAgczAgKiB4MSArIHMxICogeDIsXG4gICAgICAgIHMwICogeTEgKyBzMSAqIHkyLFxuICAgICAgICBzMCAqIHoxICsgczEgKiB6Mik7XG4gICAgfTtcbiAgfVxufTtcblxuUXVhdGVybmlvblsnWkVSTyddID0gbmV3UXVhdGVybmlvbigwLCAwLCAwLCAwKTsgLy8gVGhpcyBpcyB0aGUgYWRkaXRpdmUgaWRlbnRpdHkgUXVhdGVybmlvblxuUXVhdGVybmlvblsnT05FJ10gPSBuZXdRdWF0ZXJuaW9uKDEsIDAsIDAsIDApOyAvLyBUaGlzIGlzIHRoZSBtdWx0aXBsaWNhdGl2ZSBpZGVudGl0eSBRdWF0ZXJuaW9uXG5RdWF0ZXJuaW9uWydJJ10gPSBuZXdRdWF0ZXJuaW9uKDAsIDEsIDAsIDApO1xuUXVhdGVybmlvblsnSiddID0gbmV3UXVhdGVybmlvbigwLCAwLCAxLCAwKTtcblF1YXRlcm5pb25bJ0snXSA9IG5ld1F1YXRlcm5pb24oMCwgMCwgMCwgMSk7XG5cbi8qKlxuICogQGNvbnN0XG4gKi9cbmNvbnN0IEVQU0lMT04gPSAxZS0xNjtcblxuLyoqXG4gKiBDcmVhdGVzIHF1YXRlcm5pb24gYnkgYSByb3RhdGlvbiBnaXZlbiBhcyBheGlzLWFuZ2xlIG9yaWVudGF0aW9uXG4gKlxuICogQHBhcmFtIHtBcnJheX0gYXhpcyBUaGUgYXhpcyBhcm91bmQgd2hpY2ggdG8gcm90YXRlXG4gKiBAcGFyYW0ge251bWJlcn0gYW5nbGUgVGhlIGFuZ2xlIGluIHJhZGlhbnNcbiAqIEByZXR1cm5zIHtRdWF0ZXJuaW9ufVxuICovXG5RdWF0ZXJuaW9uWydmcm9tQXhpc0FuZ2xlJ10gPSBmdW5jdGlvbiAoYXhpcywgYW5nbGUpIHtcblxuICAvLyBRID0gW2NvcyhhbmdsZSAvIDIpLCB2ICogc2luKGFuZ2xlIC8gMildXG5cbiAgY29uc3QgYSA9IGF4aXNbMF07XG4gIGNvbnN0IGIgPSBheGlzWzFdO1xuICBjb25zdCBjID0gYXhpc1syXTtcblxuICBjb25zdCBoYWxmQW5nbGUgPSBhbmdsZSAqIDAuNTtcblxuICBjb25zdCBzaW5fMiA9IE1hdGguc2luKGhhbGZBbmdsZSk7XG4gIGNvbnN0IGNvc18yID0gTWF0aC5jb3MoaGFsZkFuZ2xlKTtcblxuICBjb25zdCBzaW5fbm9ybSA9IHNpbl8yIC8gTWF0aC5zcXJ0KGEgKiBhICsgYiAqIGIgKyBjICogYyk7XG5cbiAgcmV0dXJuIG5ld1F1YXRlcm5pb24oY29zXzIsIGEgKiBzaW5fbm9ybSwgYiAqIHNpbl9ub3JtLCBjICogc2luX25vcm0pO1xufTtcblxuLyoqXG4gKiBDYWxjdWxhdGVzIHRoZSBxdWF0ZXJuaW9uIHRvIHJvdGF0ZSB2ZWN0b3IgdSBvbnRvIHZlY3RvciB2XG4gKiBAbGluayBodHRwczovL3Jhdy5vcmcvcHJvb2YvcXVhdGVybmlvbi1mcm9tLXR3by12ZWN0b3JzL1xuICpcbiAqIEBwYXJhbSB7QXJyYXl9IHVcbiAqIEBwYXJhbSB7QXJyYXl9IHZcbiAqL1xuUXVhdGVybmlvblsnZnJvbVZlY3RvcnMnXSA9IGZ1bmN0aW9uICh1LCB2KSB7XG5cbiAgbGV0IHV4ID0gdVswXTtcbiAgbGV0IHV5ID0gdVsxXTtcbiAgbGV0IHV6ID0gdVsyXTtcblxuICBsZXQgdnggPSB2WzBdO1xuICBsZXQgdnkgPSB2WzFdO1xuICBsZXQgdnogPSB2WzJdO1xuXG4gIGNvbnN0IHVMZW4gPSBNYXRoLnNxcnQodXggKiB1eCArIHV5ICogdXkgKyB1eiAqIHV6KTtcbiAgY29uc3QgdkxlbiA9IE1hdGguc3FydCh2eCAqIHZ4ICsgdnkgKiB2eSArIHZ6ICogdnopO1xuXG4gIC8vIE5vcm1hbGl6ZSB1IGFuZCB2XG4gIGlmICh1TGVuID4gMCkgdXggLz0gdUxlbiwgdXkgLz0gdUxlbiwgdXogLz0gdUxlbjtcbiAgaWYgKHZMZW4gPiAwKSB2eCAvPSB2TGVuLCB2eSAvPSB2TGVuLCB2eiAvPSB2TGVuO1xuXG4gIC8vIENhbGN1bGF0ZSBkb3QgcHJvZHVjdCBvZiBub3JtYWxpemVkIHUgYW5kIHZcbiAgY29uc3QgZG90ID0gdXggKiB2eCArIHV5ICogdnkgKyB1eiAqIHZ6O1xuXG4gIC8vIFBhcmFsbGVsIHdoZW4gZG90ID4gMC45OTk5OTlcbiAgaWYgKGRvdCA+PSAxIC0gRVBTSUxPTikge1xuICAgIHJldHVybiBRdWF0ZXJuaW9uWydPTkUnXTtcbiAgfVxuXG4gIC8vIEFudGktUGFyYWxsZWwgKGNsb3NlIHRvIFBJKSB3aGVuIGRvdCA8IC0wLjk5OTk5OVxuICBpZiAoMSArIGRvdCA8PSBFUFNJTE9OKSB7XG5cbiAgICAvLyBSb3RhdGUgMTgwwrAgYXJvdW5kIGFueSBvcnRob2dvbmFsIHZlY3RvclxuICAgIC8vIGF4aXMgPSBsZW4oY3Jvc3MoWzEsIDAsIDBdLCB1KSkgPT0gMCA/IGNyb3NzKFswLCAxLCAwXSwgdSkgOiBjcm9zcyhbMSwgMCwgMF0sIHUpIGFuZCB0aGVyZWZvcmVcbiAgICAvLyAgICByZXR1cm4gUXVhdGVybmlvblsnZnJvbUF4aXNBbmdsZSddKE1hdGguYWJzKHV4KSA+IE1hdGguYWJzKHV6KSA/IFstdXksIHV4LCAwXSA6IFswLCAtdXosIHV5XSwgTWF0aC5QSSlcbiAgICAvLyBvciByZXR1cm4gUXVhdGVybmlvblsnZnJvbUF4aXNBbmdsZSddKE1hdGguYWJzKHV4KSA+IE1hdGguYWJzKHV6KSA/IFsgdXksLXV4LCAwXSA6IFswLCAgdXosLXV5XSwgTWF0aC5QSSlcbiAgICAvLyBvciAuLi5cblxuICAgIC8vIFNpbmNlIGZyb21BeGlzQW5nbGUoYXhpcywgUEkpID09IFF1YXRlcm5pb24oMCwgYXhpcykubm9ybWFsaXplKCksXG4gICAgaWYgKE1hdGguYWJzKHV4KSA+IE1hdGguYWJzKHV6KSkge1xuICAgICAgcmV0dXJuIG5ld05vcm1hbGl6ZWQoMCwgLXV5LCB1eCwgMCk7XG4gICAgfSBlbHNlIHtcbiAgICAgIHJldHVybiBuZXdOb3JtYWxpemVkKDAsIDAsIC11eiwgdXkpO1xuICAgIH1cbiAgfVxuXG4gIC8vIHcgPSBjcm9zcyh1LCB2KVxuICBjb25zdCB3eCA9IHV5ICogdnogLSB1eiAqIHZ5O1xuICBjb25zdCB3eSA9IHV6ICogdnggLSB1eCAqIHZ6O1xuICBjb25zdCB3eiA9IHV4ICogdnkgLSB1eSAqIHZ4O1xuXG4gIC8vIHxRfCA9IHNxcnQoKDEuMCArIGRvdCkgKiAyLjApXG4gIHJldHVybiBuZXdOb3JtYWxpemVkKDEgKyBkb3QsIHd4LCB3eSwgd3opO1xufTtcblxuLyoqXG4gKiBHZXRzIGEgc3BoZXJpY2FsIHJhbmRvbSBudW1iZXJcbiAqIEBsaW5rIGh0dHA6Ly9wbGFubmluZy5jcy51aXVjLmVkdS9ub2RlMTk4Lmh0bWxcbiAqL1xuUXVhdGVybmlvblsncmFuZG9tJ10gPSBmdW5jdGlvbiAoKSB7XG5cbiAgY29uc3QgdTEgPSBNYXRoLnJhbmRvbSgpO1xuICBjb25zdCB1MiA9IE1hdGgucmFuZG9tKCk7XG4gIGNvbnN0IHUzID0gTWF0aC5yYW5kb20oKTtcblxuICBjb25zdCBzID0gTWF0aC5zcXJ0KDEgLSB1MSk7XG4gIGNvbnN0IHQgPSBNYXRoLnNxcnQodTEpO1xuXG4gIHJldHVybiBuZXdRdWF0ZXJuaW9uKFxuICAgIHQgKiBNYXRoLmNvcygyICogTWF0aC5QSSAqIHUzKSxcbiAgICBzICogTWF0aC5zaW4oMiAqIE1hdGguUEkgKiB1MiksXG4gICAgcyAqIE1hdGguY29zKDIgKiBNYXRoLlBJICogdTIpLFxuICAgIHQgKiBNYXRoLnNpbigyICogTWF0aC5QSSAqIHUzKVxuICApO1xufTtcblxuLyoqXG4gKiBDcmVhdGVzIGEgcXVhdGVybmlvbiBieSBhIHJvdGF0aW9uIGdpdmVuIGJ5IEV1bGVyIGFuZ2xlcyAobG9naWNhbCBhcHBsaWNhdGlvbiBvcmRlciBmcm9tIGxlZnQgdG8gcmlnaHQpXG4gKlxuICogQHBhcmFtIHtudW1iZXJ9IM+IIEZpcnN0IGFuZ2xlXG4gKiBAcGFyYW0ge251bWJlcn0gzrggU2Vjb25kIGFuZ2xlXG4gKiBAcGFyYW0ge251bWJlcn0gz4YgVGhpcmQgYW5nbGVcbiAqIEBwYXJhbSB7c3RyaW5nPX0gb3JkZXIgQXhpcyBvcmRlciAoVGFpdCBCcnlhbilcbiAqIEByZXR1cm5zIHtRdWF0ZXJuaW9ufVxuICovXG5RdWF0ZXJuaW9uWydmcm9tRXVsZXJMb2dpY2FsJ10gPSBmdW5jdGlvbiAoz4gsIM64LCDPhiwgb3JkZXIpIHtcblxuICByZXR1cm4gUXVhdGVybmlvblsnZnJvbUV1bGVyJ10oz4YsIM64LCDPiCwgb3JkZXIgIT09IHVuZGVmaW5lZCA/IG9yZGVyWzJdICsgb3JkZXJbMV0gKyBvcmRlclswXSA6IG9yZGVyKTtcbn07XG5cbi8qKlxuICogQ3JlYXRlcyBhIHF1YXRlcm5pb24gYnkgYSByb3RhdGlvbiBnaXZlbiBieSBFdWxlciBhbmdsZXMgKG11bHRpcGxpY2F0aW9uIG9yZGVyIGZyb20gcmlnaHQgdG8gbGVmdClcbiAqXG4gKiBAcGFyYW0ge251bWJlcn0gz4YgRmlyc3QgYW5nbGVcbiAqIEBwYXJhbSB7bnVtYmVyfSDOuCBTZWNvbmQgYW5nbGVcbiAqIEBwYXJhbSB7bnVtYmVyfSDPiCBUaGlyZCBhbmdsZVxuICogQHBhcmFtIHtzdHJpbmc9fSBvcmRlciBBeGlzIG9yZGVyIChUYWl0IEJyeWFuKVxuICogQHJldHVybnMge1F1YXRlcm5pb259XG4gKi9cblF1YXRlcm5pb25bJ2Zyb21FdWxlciddID0gZnVuY3Rpb24gKM+GLCDOuCwgz4gsIG9yZGVyKSB7XG5cbiAgY29uc3QgX3ggPSDPhiAqIDAuNTtcbiAgY29uc3QgX3kgPSDOuCAqIDAuNTtcbiAgY29uc3QgX3ogPSDPiCAqIDAuNTtcblxuICBjb25zdCBjWCA9IE1hdGguY29zKF94KTtcbiAgY29uc3QgY1kgPSBNYXRoLmNvcyhfeSk7XG4gIGNvbnN0IGNaID0gTWF0aC5jb3MoX3opO1xuXG4gIGNvbnN0IHNYID0gTWF0aC5zaW4oX3gpO1xuICBjb25zdCBzWSA9IE1hdGguc2luKF95KTtcbiAgY29uc3Qgc1ogPSBNYXRoLnNpbihfeik7XG5cbiAgaWYgKG9yZGVyID09PSB1bmRlZmluZWQgfHwgb3JkZXIgPT09ICdaWFknKSB7XG4gICAgLy8gYXhpc0FuZ2xlKFswLCAwLCAxXSwgz4YpICogYXhpc0FuZ2xlKFsxLCAwLCAwXSwgzrgpICogYXhpc0FuZ2xlKFswLCAxLCAwXSwgz4gpXG4gICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oXG4gICAgICBjWCAqIGNZICogY1ogLSBzWCAqIHNZICogc1osXG4gICAgICBzWSAqIGNYICogY1ogLSBzWCAqIHNaICogY1ksXG4gICAgICBzWCAqIHNZICogY1ogKyBzWiAqIGNYICogY1ksXG4gICAgICBzWCAqIGNZICogY1ogKyBzWSAqIHNaICogY1gpO1xuICB9XG5cbiAgaWYgKG9yZGVyID09PSAnWFlaJyB8fCBvcmRlciA9PT0gJ1JQWScpIHsgLy8gcm9sbCBhcm91bmQgWCwgcGl0Y2ggYXJvdW5kIFksIHlhdyBhcm91bmQgWlxuICAgIC8vIGF4aXNBbmdsZShbMSwgMCwgMF0sIM+GKSAqIGF4aXNBbmdsZShbMCwgMSwgMF0sIM64KSAqIGF4aXNBbmdsZShbMCwgMCwgMV0sIM+IKVxuICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKFxuICAgICAgY1ggKiBjWSAqIGNaIC0gc1ggKiBzWSAqIHNaLFxuICAgICAgc1ggKiBjWSAqIGNaICsgc1kgKiBzWiAqIGNYLFxuICAgICAgc1kgKiBjWCAqIGNaIC0gc1ggKiBzWiAqIGNZLFxuICAgICAgc1ggKiBzWSAqIGNaICsgc1ogKiBjWCAqIGNZKTtcbiAgfVxuXG4gIGlmIChvcmRlciA9PT0gJ1lYWicpIHsgLy8gZGV2aWNlb3JpZW50YXRpb25cbiAgICAvLyBheGlzQW5nbGUoWzAsIDEsIDBdLCDPhikgKiBheGlzQW5nbGUoWzEsIDAsIDBdLCDOuCkgKiBheGlzQW5nbGUoWzAsIDAsIDFdLCDPiClcbiAgICByZXR1cm4gbmV3UXVhdGVybmlvbihcbiAgICAgIHNYICogc1kgKiBzWiArIGNYICogY1kgKiBjWixcbiAgICAgIHNYICogc1ogKiBjWSArIHNZICogY1ggKiBjWixcbiAgICAgIHNYICogY1kgKiBjWiAtIHNZICogc1ogKiBjWCxcbiAgICAgIHNaICogY1ggKiBjWSAtIHNYICogc1kgKiBjWik7XG4gIH1cblxuICBpZiAob3JkZXIgPT09ICdaWVgnIHx8IG9yZGVyID09PSAnWVBSJykge1xuICAgIC8vIGF4aXNBbmdsZShbMCwgMCwgMV0sIM+GKSAqIGF4aXNBbmdsZShbMCwgMSwgMF0sIM64KSAqIGF4aXNBbmdsZShbMSwgMCwgMF0sIM+IKVxuICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKFxuICAgICAgc1ggKiBzWSAqIHNaICsgY1ggKiBjWSAqIGNaLFxuICAgICAgc1ogKiBjWCAqIGNZIC0gc1ggKiBzWSAqIGNaLFxuICAgICAgc1ggKiBzWiAqIGNZICsgc1kgKiBjWCAqIGNaLFxuICAgICAgc1ggKiBjWSAqIGNaIC0gc1kgKiBzWiAqIGNYKTtcbiAgfVxuXG4gIGlmIChvcmRlciA9PT0gJ1laWCcpIHtcbiAgICAvLyBheGlzQW5nbGUoWzAsIDEsIDBdLCDPhikgKiBheGlzQW5nbGUoWzAsIDAsIDFdLCDOuCkgKiBheGlzQW5nbGUoWzEsIDAsIDBdLCDPiClcbiAgICByZXR1cm4gbmV3UXVhdGVybmlvbihcbiAgICAgIGNYICogY1kgKiBjWiAtIHNYICogc1kgKiBzWixcbiAgICAgIHNYICogc1kgKiBjWiArIHNaICogY1ggKiBjWSxcbiAgICAgIHNYICogY1kgKiBjWiArIHNZICogc1ogKiBjWCxcbiAgICAgIHNZICogY1ggKiBjWiAtIHNYICogc1ogKiBjWSk7XG4gIH1cblxuICBpZiAob3JkZXIgPT09ICdYWlknKSB7XG4gICAgLy8gYXhpc0FuZ2xlKFsxLCAwLCAwXSwgz4YpICogYXhpc0FuZ2xlKFswLCAwLCAxXSwgzrgpICogYXhpc0FuZ2xlKFswLCAxLCAwXSwgz4gpXG4gICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oXG4gICAgICBzWCAqIHNZICogc1ogKyBjWCAqIGNZICogY1osXG4gICAgICBzWCAqIGNZICogY1ogLSBzWSAqIHNaICogY1gsXG4gICAgICBzWiAqIGNYICogY1kgLSBzWCAqIHNZICogY1osXG4gICAgICBzWCAqIHNaICogY1kgKyBzWSAqIGNYICogY1opO1xuICB9XG5cbiAgaWYgKG9yZGVyID09PSAnWllaJykge1xuICAgIC8vIGF4aXNBbmdsZShbMCwgMCwgMV0sIM+GKSAqIGF4aXNBbmdsZShbMCwgMSwgMF0sIM64KSAqIGF4aXNBbmdsZShbMCwgMCwgMV0sIM+IKVxuICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKFxuICAgICAgY1ggKiBjWSAqIGNaIC0gc1ggKiBzWiAqIGNZLFxuICAgICAgc1kgKiBzWiAqIGNYIC0gc1ggKiBzWSAqIGNaLFxuICAgICAgc1ggKiBzWSAqIHNaICsgc1kgKiBjWCAqIGNaLFxuICAgICAgc1ggKiBjWSAqIGNaICsgc1ogKiBjWCAqIGNZKTtcbiAgfVxuXG4gIGlmIChvcmRlciA9PT0gJ1pYWicpIHtcbiAgICAvLyBheGlzQW5nbGUoWzAsIDAsIDFdLCDPhikgKiBheGlzQW5nbGUoWzEsIDAsIDBdLCDOuCkgKiBheGlzQW5nbGUoWzAsIDAsIDFdLCDPiClcbiAgICByZXR1cm4gbmV3UXVhdGVybmlvbihcbiAgICAgIGNYICogY1kgKiBjWiAtIHNYICogc1ogKiBjWSxcbiAgICAgIHNYICogc1kgKiBzWiArIHNZICogY1ggKiBjWixcbiAgICAgIHNYICogc1kgKiBjWiAtIHNZICogc1ogKiBjWCxcbiAgICAgIHNYICogY1kgKiBjWiArIHNaICogY1ggKiBjWSk7XG4gIH1cblxuICBpZiAob3JkZXIgPT09ICdZWFknKSB7XG4gICAgLy8gYXhpc0FuZ2xlKFswLCAxLCAwXSwgz4YpICogYXhpc0FuZ2xlKFsxLCAwLCAwXSwgzrgpICogYXhpc0FuZ2xlKFswLCAxLCAwXSwgz4gpXG4gICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oXG4gICAgICBjWCAqIGNZICogY1ogLSBzWCAqIHNaICogY1ksXG4gICAgICBzWCAqIHNZICogc1ogKyBzWSAqIGNYICogY1osXG4gICAgICBzWCAqIGNZICogY1ogKyBzWiAqIGNYICogY1ksXG4gICAgICBzWSAqIHNaICogY1ggLSBzWCAqIHNZICogY1opO1xuICB9XG5cbiAgaWYgKG9yZGVyID09PSAnWVpZJykge1xuICAgIC8vIGF4aXNBbmdsZShbMCwgMSwgMF0sIM+GKSAqIGF4aXNBbmdsZShbMCwgMCwgMV0sIM64KSAqIGF4aXNBbmdsZShbMCwgMSwgMF0sIM+IKVxuICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKFxuICAgICAgY1ggKiBjWSAqIGNaIC0gc1ggKiBzWiAqIGNZLFxuICAgICAgc1ggKiBzWSAqIGNaIC0gc1kgKiBzWiAqIGNYLFxuICAgICAgc1ggKiBjWSAqIGNaICsgc1ogKiBjWCAqIGNZLFxuICAgICAgc1ggKiBzWSAqIHNaICsgc1kgKiBjWCAqIGNaKTtcbiAgfVxuXG4gIGlmIChvcmRlciA9PT0gJ1hZWCcpIHtcbiAgICAvLyBheGlzQW5nbGUoWzEsIDAsIDBdLCDPhikgKiBheGlzQW5nbGUoWzAsIDEsIDBdLCDOuCkgKiBheGlzQW5nbGUoWzEsIDAsIDBdLCDPiClcbiAgICByZXR1cm4gbmV3UXVhdGVybmlvbihcbiAgICAgIGNYICogY1kgKiBjWiAtIHNYICogc1ogKiBjWSxcbiAgICAgIHNYICogY1kgKiBjWiArIHNaICogY1ggKiBjWSxcbiAgICAgIHNYICogc1kgKiBzWiArIHNZICogY1ggKiBjWixcbiAgICAgIHNYICogc1kgKiBjWiAtIHNZICogc1ogKiBjWCk7XG4gIH1cblxuICBpZiAob3JkZXIgPT09ICdYWlgnKSB7XG4gICAgLy8gYXhpc0FuZ2xlKFsxLCAwLCAwXSwgz4YpICogYXhpc0FuZ2xlKFswLCAwLCAxXSwgzrgpICogYXhpc0FuZ2xlKFsxLCAwLCAwXSwgz4gpXG4gICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oXG4gICAgICBjWCAqIGNZICogY1ogLSBzWCAqIHNaICogY1ksXG4gICAgICBzWCAqIGNZICogY1ogKyBzWiAqIGNYICogY1ksXG4gICAgICBzWSAqIHNaICogY1ggLSBzWCAqIHNZICogY1osXG4gICAgICBzWCAqIHNZICogc1ogKyBzWSAqIGNYICogY1opO1xuICB9XG4gIHJldHVybiBudWxsO1xufTtcblxuLyoqXG4gKiBDcmVhdGVzIGEgcXVhdGVybmlvbiBieSBhIHJvdGF0aW9uIG1hdHJpeFxuICpcbiAqIEBwYXJhbSB7QXJyYXl9IG1hdHJpeFxuICogQHJldHVybnMge1F1YXRlcm5pb259XG4gKi9cblF1YXRlcm5pb25bJ2Zyb21NYXRyaXgnXSA9IGZ1bmN0aW9uIChtYXRyaXgpIHtcblxuICBsZXQgbTAwLCBtMDEsIG0wMiwgbTEwLCBtMTEsIG0xMiwgbTIwLCBtMjEsIG0yMjtcblxuICBpZiAobWF0cml4Lmxlbmd0aCA9PT0gOSkge1xuICAgIG0wMCA9IG1hdHJpeFswXTtcbiAgICBtMDEgPSBtYXRyaXhbMV07XG4gICAgbTAyID0gbWF0cml4WzJdO1xuXG4gICAgbTEwID0gbWF0cml4WzNdO1xuICAgIG0xMSA9IG1hdHJpeFs0XTtcbiAgICBtMTIgPSBtYXRyaXhbNV07XG5cbiAgICBtMjAgPSBtYXRyaXhbNl07XG4gICAgbTIxID0gbWF0cml4WzddO1xuICAgIG0yMiA9IG1hdHJpeFs4XTtcblxuICB9IGVsc2Uge1xuICAgIG0wMCA9IG1hdHJpeFswXVswXTtcbiAgICBtMDEgPSBtYXRyaXhbMF1bMV07XG4gICAgbTAyID0gbWF0cml4WzBdWzJdO1xuXG4gICAgbTEwID0gbWF0cml4WzFdWzBdO1xuICAgIG0xMSA9IG1hdHJpeFsxXVsxXTtcbiAgICBtMTIgPSBtYXRyaXhbMV1bMl07XG5cbiAgICBtMjAgPSBtYXRyaXhbMl1bMF07XG4gICAgbTIxID0gbWF0cml4WzJdWzFdO1xuICAgIG0yMiA9IG1hdHJpeFsyXVsyXTtcbiAgfVxuXG4gIGNvbnN0IHRyID0gbTAwICsgbTExICsgbTIyOyAvLyAyICogdyA9IHNxcnQoMSArIHRyKVxuXG4gIC8vIENob29zZSB0aGUgZWxlbWVudCB3aXRoIHRoZSBiaWdnZXN0IHZhbHVlIG9uIHRoZSBkaWFnb25hbFxuXG4gIGlmICh0ciA+IDApIHsgLy8gaWYgdHJhY2UgaXMgcG9zaXRpdmUgdGhlbiBcIndcIiBpcyBiaWdnZXN0IGNvbXBvbmVudFxuICAgIC8vIHxRfCA9IDIgKiBzcXJ0KDEgKyB0cikgPSA0d1xuICAgIHJldHVybiBuZXdOb3JtYWxpemVkKFxuICAgICAgdHIgKyAxLjAsXG4gICAgICBtMjEgLSBtMTIsXG4gICAgICBtMDIgLSBtMjAsXG4gICAgICBtMTAgLSBtMDEpO1xuICB9IGVsc2UgaWYgKG0wMCA+IG0xMSAmJiBtMDAgPiBtMjIpIHtcbiAgICAvLyB8UXwgPSAyICogc3FydCgxLjAgKyBtMDAgLSBtMTEgLSBtMjIpID0gNHhcbiAgICByZXR1cm4gbmV3Tm9ybWFsaXplZChcbiAgICAgIG0yMSAtIG0xMixcbiAgICAgIDEuMCArIG0wMCAtIG0xMSAtIG0yMixcbiAgICAgIG0wMSArIG0xMCxcbiAgICAgIG0wMiArIG0yMCk7XG4gIH0gZWxzZSBpZiAobTExID4gbTIyKSB7XG4gICAgLy8gfFF8ID0gMiAqIHNxcnQoMS4wICsgbTExIC0gbTAwIC0gbTIyKSA9IDR5XG4gICAgcmV0dXJuIG5ld05vcm1hbGl6ZWQoXG4gICAgICBtMDIgLSBtMjAsXG4gICAgICBtMDEgKyBtMTAsXG4gICAgICAxLjAgKyBtMTEgLSBtMDAgLSBtMjIsXG4gICAgICBtMTIgKyBtMjEpO1xuICB9IGVsc2Uge1xuICAgIC8vIHxRfCA9IDIgKiBzcXJ0KDEuMCArIG0yMiAtIG0wMCAtIG0xMSkgPSA0elxuICAgIHJldHVybiBuZXdOb3JtYWxpemVkKFxuICAgICAgbTEwIC0gbTAxLFxuICAgICAgbTAyICsgbTIwLFxuICAgICAgbTEyICsgbTIxLFxuICAgICAgMS4wICsgbTIyIC0gbTAwIC0gbTExKTtcbiAgfVxufTtcblxuT2JqZWN0LmRlZmluZVByb3BlcnR5KFF1YXRlcm5pb24sIFwiX19lc01vZHVsZVwiLCB7ICd2YWx1ZSc6IHRydWUgfSk7XG5RdWF0ZXJuaW9uWydkZWZhdWx0J10gPSBRdWF0ZXJuaW9uO1xuUXVhdGVybmlvblsnUXVhdGVybmlvbiddID0gUXVhdGVybmlvbjtcbm1vZHVsZVsnZXhwb3J0cyddID0gUXVhdGVybmlvbjtcbiIsIlwidXNlIHN0cmljdFwiO1xuT2JqZWN0LmRlZmluZVByb3BlcnR5KGV4cG9ydHMsIFwiX19lc01vZHVsZVwiLCB7IHZhbHVlOiB0cnVlIH0pO1xuZXhwb3J0cy5hbGlhc2VzID0gdm9pZCAwO1xuZXhwb3J0cy5hbGlhc2VzID0ge1xuICAgIFwibVBlbHZpc1wiOiBcImhpcFwiLFxuICAgIFwibVNwaW5lMVwiOiBcIm1TcGluZTFcIixcbiAgICBcIm1TcGluZTJcIjogXCJtU3BpbmUyXCIsXG4gICAgXCJtVG9yc29cIjogXCJhYmRvbWVuXCIsXG4gICAgXCJtU3BpbmUzXCI6IFwibVNwaW5lM1wiLFxuICAgIFwibVNwaW5lNFwiOiBcIm1TcGluZTRcIixcbiAgICBcIm1DaGVzdFwiOiBcImNoZXN0XCIsXG4gICAgXCJtTmVja1wiOiBcIm5lY2tcIixcbiAgICBcIm1IZWFkXCI6IFwiaGVhZFwiLFxuICAgIFwibVNrdWxsXCI6IFwiZmlndXJlSGFpclwiLFxuICAgIFwibUV5ZVJpZ2h0XCI6IFwibUV5ZVJpZ2h0XCIsXG4gICAgXCJtRXllTGVmdFwiOiBcIm1FeWVMZWZ0XCIsXG4gICAgXCJtRmFjZVJvb3RcIjogXCJtRmFjZVJvb3RcIixcbiAgICBcIm1GYWNlRXllQWx0UmlnaHRcIjogXCJtRmFjZUV5ZUFsdFJpZ2h0XCIsXG4gICAgXCJtRmFjZUV5ZUFsdExlZnRcIjogXCJtRmFjZUV5ZUFsdExlZnRcIixcbiAgICBcIm1GYWNlRm9yZWhlYWRMZWZ0XCI6IFwibUZhY2VGb3JlaGVhZExlZnRcIixcbiAgICBcIm1GYWNlRm9yZWhlYWRSaWdodFwiOiBcIm1GYWNlRm9yZWhlYWRSaWdodFwiLFxuICAgIFwibUZhY2VFeWVicm93T3V0ZXJMZWZ0XCI6IFwibUZhY2VFeWVicm93T3V0ZXJMZWZ0XCIsXG4gICAgXCJtRmFjZUV5ZWJyb3dDZW50ZXJMZWZ0XCI6IFwibUZhY2VFeWVicm93Q2VudGVyTGVmdFwiLFxuICAgIFwibUZhY2VFeWVicm93SW5uZXJMZWZ0XCI6IFwibUZhY2VFeWVicm93SW5uZXJMZWZ0XCIsXG4gICAgXCJtRmFjZUV5ZWJyb3dPdXRlclJpZ2h0XCI6IFwibUZhY2VFeWVicm93T3V0ZXJSaWdodFwiLFxuICAgIFwibUZhY2VFeWVicm93Q2VudGVyUmlnaHRcIjogXCJtRmFjZUV5ZWJyb3dDZW50ZXJSaWdodFwiLFxuICAgIFwibUZhY2VFeWVicm93SW5uZXJSaWdodFwiOiBcIm1GYWNlRXllYnJvd0lubmVyUmlnaHRcIixcbiAgICBcIm1GYWNlRXllTGlkVXBwZXJMZWZ0XCI6IFwibUZhY2VFeWVMaWRVcHBlckxlZnRcIixcbiAgICBcIm1GYWNlRXllTGlkTG93ZXJMZWZ0XCI6IFwibUZhY2VFeWVMaWRMb3dlckxlZnRcIixcbiAgICBcIm1GYWNlRXllTGlkVXBwZXJSaWdodFwiOiBcIm1GYWNlRXllTGlkVXBwZXJSaWdodFwiLFxuICAgIFwibUZhY2VFeWVMaWRMb3dlclJpZ2h0XCI6IFwibUZhY2VFeWVMaWRMb3dlclJpZ2h0XCIsXG4gICAgXCJtRmFjZUVhcjFMZWZ0XCI6IFwibUZhY2VFYXIxTGVmdFwiLFxuICAgIFwibUZhY2VFYXIyTGVmdFwiOiBcIm1GYWNlRWFyMkxlZnRcIixcbiAgICBcIm1GYWNlRWFyMVJpZ2h0XCI6IFwibUZhY2VFYXIxUmlnaHRcIixcbiAgICBcIm1GYWNlRWFyMlJpZ2h0XCI6IFwibUZhY2VFYXIyUmlnaHRcIixcbiAgICBcIm1GYWNlTm9zZUxlZnRcIjogXCJtRmFjZU5vc2VMZWZ0XCIsXG4gICAgXCJtRmFjZU5vc2VDZW50ZXJcIjogXCJtRmFjZU5vc2VDZW50ZXJcIixcbiAgICBcIm1GYWNlTm9zZVJpZ2h0XCI6IFwibUZhY2VOb3NlUmlnaHRcIixcbiAgICBcIm1GYWNlQ2hlZWtMb3dlckxlZnRcIjogXCJtRmFjZUNoZWVrTG93ZXJMZWZ0XCIsXG4gICAgXCJtRmFjZUNoZWVrVXBwZXJMZWZ0XCI6IFwibUZhY2VDaGVla1VwcGVyTGVmdFwiLFxuICAgIFwibUZhY2VDaGVla0xvd2VyUmlnaHRcIjogXCJtRmFjZUNoZWVrTG93ZXJSaWdodFwiLFxuICAgIFwibUZhY2VDaGVla1VwcGVyUmlnaHRcIjogXCJtRmFjZUNoZWVrVXBwZXJSaWdodFwiLFxuICAgIFwibUZhY2VKYXdcIjogXCJtRmFjZUphd1wiLFxuICAgIFwibUZhY2VDaGluXCI6IFwibUZhY2VDaGluXCIsXG4gICAgXCJtRmFjZVRlZXRoTG93ZXJcIjogXCJtRmFjZVRlZXRoTG93ZXJcIixcbiAgICBcIm1GYWNlTGlwTG93ZXJMZWZ0XCI6IFwibUZhY2VMaXBMb3dlckxlZnRcIixcbiAgICBcIm1GYWNlTGlwTG93ZXJSaWdodFwiOiBcIm1GYWNlTGlwTG93ZXJSaWdodFwiLFxuICAgIFwibUZhY2VMaXBMb3dlckNlbnRlclwiOiBcIm1GYWNlTGlwTG93ZXJDZW50ZXJcIixcbiAgICBcIm1GYWNlVG9uZ3VlQmFzZVwiOiBcIm1GYWNlVG9uZ3VlQmFzZVwiLFxuICAgIFwibUZhY2VUb25ndWVUaXBcIjogXCJtRmFjZVRvbmd1ZVRpcFwiLFxuICAgIFwibUZhY2VKYXdTaGFwZXJcIjogXCJtRmFjZUphd1NoYXBlclwiLFxuICAgIFwibUZhY2VGb3JlaGVhZENlbnRlclwiOiBcIm1GYWNlRm9yZWhlYWRDZW50ZXJcIixcbiAgICBcIm1GYWNlTm9zZUJhc2VcIjogXCJtRmFjZU5vc2VCYXNlXCIsXG4gICAgXCJtRmFjZVRlZXRoVXBwZXJcIjogXCJtRmFjZVRlZXRoVXBwZXJcIixcbiAgICBcIm1GYWNlTGlwVXBwZXJMZWZ0XCI6IFwibUZhY2VMaXBVcHBlckxlZnRcIixcbiAgICBcIm1GYWNlTGlwVXBwZXJSaWdodFwiOiBcIm1GYWNlTGlwVXBwZXJSaWdodFwiLFxuICAgIFwibUZhY2VMaXBDb3JuZXJMZWZ0XCI6IFwibUZhY2VMaXBDb3JuZXJMZWZ0XCIsXG4gICAgXCJtRmFjZUxpcENvcm5lclJpZ2h0XCI6IFwibUZhY2VMaXBDb3JuZXJSaWdodFwiLFxuICAgIFwibUZhY2VMaXBVcHBlckNlbnRlclwiOiBcIm1GYWNlTGlwVXBwZXJDZW50ZXJcIixcbiAgICBcIm1GYWNlRXllY29ybmVySW5uZXJMZWZ0XCI6IFwibUZhY2VFeWVjb3JuZXJJbm5lckxlZnRcIixcbiAgICBcIm1GYWNlRXllY29ybmVySW5uZXJSaWdodFwiOiBcIm1GYWNlRXllY29ybmVySW5uZXJSaWdodFwiLFxuICAgIFwibUZhY2VOb3NlQnJpZGdlXCI6IFwibUZhY2VOb3NlQnJpZGdlXCIsXG4gICAgXCJtQ29sbGFyTGVmdFwiOiBcImxDb2xsYXJcIixcbiAgICBcIm1TaG91bGRlckxlZnRcIjogXCJsU2hsZHJcIixcbiAgICBcIm1FbGJvd0xlZnRcIjogXCJsRm9yZUFybVwiLFxuICAgIFwibVdyaXN0TGVmdFwiOiBcImxIYW5kXCIsXG4gICAgXCJtSGFuZE1pZGRsZTFMZWZ0XCI6IFwibUhhbmRNaWRkbGUxTGVmdFwiLFxuICAgIFwibUhhbmRNaWRkbGUyTGVmdFwiOiBcIm1IYW5kTWlkZGxlMkxlZnRcIixcbiAgICBcIm1IYW5kTWlkZGxlM0xlZnRcIjogXCJtSGFuZE1pZGRsZTNMZWZ0XCIsXG4gICAgXCJtSGFuZEluZGV4MUxlZnRcIjogXCJtSGFuZEluZGV4MUxlZnRcIixcbiAgICBcIm1IYW5kSW5kZXgyTGVmdFwiOiBcIm1IYW5kSW5kZXgyTGVmdFwiLFxuICAgIFwibUhhbmRJbmRleDNMZWZ0XCI6IFwibUhhbmRJbmRleDNMZWZ0XCIsXG4gICAgXCJtSGFuZFJpbmcxTGVmdFwiOiBcIm1IYW5kUmluZzFMZWZ0XCIsXG4gICAgXCJtSGFuZFJpbmcyTGVmdFwiOiBcIm1IYW5kUmluZzJMZWZ0XCIsXG4gICAgXCJtSGFuZFJpbmczTGVmdFwiOiBcIm1IYW5kUmluZzNMZWZ0XCIsXG4gICAgXCJtSGFuZFBpbmt5MUxlZnRcIjogXCJtSGFuZFBpbmt5MUxlZnRcIixcbiAgICBcIm1IYW5kUGlua3kyTGVmdFwiOiBcIm1IYW5kUGlua3kyTGVmdFwiLFxuICAgIFwibUhhbmRQaW5reTNMZWZ0XCI6IFwibUhhbmRQaW5reTNMZWZ0XCIsXG4gICAgXCJtSGFuZFRodW1iMUxlZnRcIjogXCJtSGFuZFRodW1iMUxlZnRcIixcbiAgICBcIm1IYW5kVGh1bWIyTGVmdFwiOiBcIm1IYW5kVGh1bWIyTGVmdFwiLFxuICAgIFwibUhhbmRUaHVtYjNMZWZ0XCI6IFwibUhhbmRUaHVtYjNMZWZ0XCIsXG4gICAgXCJtQ29sbGFyUmlnaHRcIjogXCJyQ29sbGFyXCIsXG4gICAgXCJtU2hvdWxkZXJSaWdodFwiOiBcInJTaGxkclwiLFxuICAgIFwibUVsYm93UmlnaHRcIjogXCJyRm9yZUFybVwiLFxuICAgIFwibVdyaXN0UmlnaHRcIjogXCJySGFuZFwiLFxuICAgIFwibUhhbmRNaWRkbGUxUmlnaHRcIjogXCJtSGFuZE1pZGRsZTFSaWdodFwiLFxuICAgIFwibUhhbmRNaWRkbGUyUmlnaHRcIjogXCJtSGFuZE1pZGRsZTJSaWdodFwiLFxuICAgIFwibUhhbmRNaWRkbGUzUmlnaHRcIjogXCJtSGFuZE1pZGRsZTNSaWdodFwiLFxuICAgIFwibUhhbmRJbmRleDFSaWdodFwiOiBcIm1IYW5kSW5kZXgxUmlnaHRcIixcbiAgICBcIm1IYW5kSW5kZXgyUmlnaHRcIjogXCJtSGFuZEluZGV4MlJpZ2h0XCIsXG4gICAgXCJtSGFuZEluZGV4M1JpZ2h0XCI6IFwibUhhbmRJbmRleDNSaWdodFwiLFxuICAgIFwibUhhbmRSaW5nMVJpZ2h0XCI6IFwibUhhbmRSaW5nMVJpZ2h0XCIsXG4gICAgXCJtSGFuZFJpbmcyUmlnaHRcIjogXCJtSGFuZFJpbmcyUmlnaHRcIixcbiAgICBcIm1IYW5kUmluZzNSaWdodFwiOiBcIm1IYW5kUmluZzNSaWdodFwiLFxuICAgIFwibUhhbmRQaW5reTFSaWdodFwiOiBcIm1IYW5kUGlua3kxUmlnaHRcIixcbiAgICBcIm1IYW5kUGlua3kyUmlnaHRcIjogXCJtSGFuZFBpbmt5MlJpZ2h0XCIsXG4gICAgXCJtSGFuZFBpbmt5M1JpZ2h0XCI6IFwibUhhbmRQaW5reTNSaWdodFwiLFxuICAgIFwibUhhbmRUaHVtYjFSaWdodFwiOiBcIm1IYW5kVGh1bWIxUmlnaHRcIixcbiAgICBcIm1IYW5kVGh1bWIyUmlnaHRcIjogXCJtSGFuZFRodW1iMlJpZ2h0XCIsXG4gICAgXCJtSGFuZFRodW1iM1JpZ2h0XCI6IFwibUhhbmRUaHVtYjNSaWdodFwiLFxuICAgIFwibVdpbmdzUm9vdFwiOiBcIm1XaW5nc1Jvb3RcIixcbiAgICBcIm1XaW5nMUxlZnRcIjogXCJtV2luZzFMZWZ0XCIsXG4gICAgXCJtV2luZzJMZWZ0XCI6IFwibVdpbmcyTGVmdFwiLFxuICAgIFwibVdpbmczTGVmdFwiOiBcIm1XaW5nM0xlZnRcIixcbiAgICBcIm1XaW5nNExlZnRcIjogXCJtV2luZzRMZWZ0XCIsXG4gICAgXCJtV2luZzRGYW5MZWZ0XCI6IFwibVdpbmc0RmFuTGVmdFwiLFxuICAgIFwibVdpbmcxUmlnaHRcIjogXCJtV2luZzFSaWdodFwiLFxuICAgIFwibVdpbmcyUmlnaHRcIjogXCJtV2luZzJSaWdodFwiLFxuICAgIFwibVdpbmczUmlnaHRcIjogXCJtV2luZzNSaWdodFwiLFxuICAgIFwibVdpbmc0UmlnaHRcIjogXCJtV2luZzRSaWdodFwiLFxuICAgIFwibVdpbmc0RmFuUmlnaHRcIjogXCJtV2luZzRGYW5SaWdodFwiLFxuICAgIFwibUhpcFJpZ2h0XCI6IFwiclRoaWdoXCIsXG4gICAgXCJtS25lZVJpZ2h0XCI6IFwiclNoaW5cIixcbiAgICBcIm1BbmtsZVJpZ2h0XCI6IFwickZvb3RcIixcbiAgICBcIm1Gb290UmlnaHRcIjogXCJtRm9vdFJpZ2h0XCIsXG4gICAgXCJtVG9lUmlnaHRcIjogXCJtVG9lUmlnaHRcIixcbiAgICBcIm1IaXBMZWZ0XCI6IFwibFRoaWdoXCIsXG4gICAgXCJtS25lZUxlZnRcIjogXCJsU2hpblwiLFxuICAgIFwibUFua2xlTGVmdFwiOiBcImxGb290XCIsXG4gICAgXCJtRm9vdExlZnRcIjogXCJtRm9vdExlZnRcIixcbiAgICBcIm1Ub2VMZWZ0XCI6IFwibVRvZUxlZnRcIixcbiAgICBcIm1UYWlsMVwiOiBcIm1UYWlsMVwiLFxuICAgIFwibVRhaWwyXCI6IFwibVRhaWwyXCIsXG4gICAgXCJtVGFpbDNcIjogXCJtVGFpbDNcIixcbiAgICBcIm1UYWlsNFwiOiBcIm1UYWlsNFwiLFxuICAgIFwibVRhaWw1XCI6IFwibVRhaWw1XCIsXG4gICAgXCJtVGFpbDZcIjogXCJtVGFpbDZcIixcbiAgICBcIm1Hcm9pblwiOiBcIm1Hcm9pblwiLFxuICAgIFwibUhpbmRMaW1ic1Jvb3RcIjogXCJtSGluZExpbWJzUm9vdFwiLFxuICAgIFwibUhpbmRMaW1iMUxlZnRcIjogXCJtSGluZExpbWIxTGVmdFwiLFxuICAgIFwibUhpbmRMaW1iMkxlZnRcIjogXCJtSGluZExpbWIyTGVmdFwiLFxuICAgIFwibUhpbmRMaW1iM0xlZnRcIjogXCJtSGluZExpbWIzTGVmdFwiLFxuICAgIFwibUhpbmRMaW1iNExlZnRcIjogXCJtSGluZExpbWI0TGVmdFwiLFxuICAgIFwibUhpbmRMaW1iMVJpZ2h0XCI6IFwibUhpbmRMaW1iMVJpZ2h0XCIsXG4gICAgXCJtSGluZExpbWIyUmlnaHRcIjogXCJtSGluZExpbWIyUmlnaHRcIixcbiAgICBcIm1IaW5kTGltYjNSaWdodFwiOiBcIm1IaW5kTGltYjNSaWdodFwiLFxuICAgIFwibUhpbmRMaW1iNFJpZ2h0XCI6IFwibUhpbmRMaW1iNFJpZ2h0XCJcbn07XG4iLCJcInVzZSBzdHJpY3RcIjtcbk9iamVjdC5kZWZpbmVQcm9wZXJ0eShleHBvcnRzLCBcIl9fZXNNb2R1bGVcIiwgeyB2YWx1ZTogdHJ1ZSB9KTtcbmV4cG9ydHMudmlzaXROb2RlID0gdmlzaXROb2RlO1xuZXhwb3J0cy5zZXJpYWxpemVCVkggPSBzZXJpYWxpemVCVkg7XG5leHBvcnRzLnRvQlZIID0gdG9CVkg7XG5jb25zdCB1dGlsc18xID0gcmVxdWlyZShcIi4vdXRpbHNcIik7XG5jb25zdCBoaWVyYXJjaHlfMSA9IHJlcXVpcmUoXCIuL2hpZXJhcmNoeVwiKTtcbmNvbnN0IGFsaWFzZXNfMSA9IHJlcXVpcmUoXCIuL2FsaWFzZXNcIik7XG5jb25zdCBxdWF0ZXJuaW9uXzEgPSByZXF1aXJlKFwicXVhdGVybmlvblwiKTtcbmZ1bmN0aW9uIG9mZnNldFRvU3RyaW5nKG9mZnNldCwgZGlnaXRzKSB7XG4gICAgcmV0dXJuIFwiT0ZGU0VUIFwiICsgKDAsIHV0aWxzXzEuZmxvYXRUb1N0cmluZykob2Zmc2V0LngsIGRpZ2l0cykgKyBcIiBcIiArICgwLCB1dGlsc18xLmZsb2F0VG9TdHJpbmcpKG9mZnNldC55LCBkaWdpdHMpICsgXCIgXCIgKyAoMCwgdXRpbHNfMS5mbG9hdFRvU3RyaW5nKShvZmZzZXQueiwgZGlnaXRzKTtcbn1cbmZ1bmN0aW9uIGNoYW5uZWxzU3RyaW5nKG5vZGUpIHtcbiAgICBpZiAoIW5vZGUuY2hhbm5lbHMpIHtcbiAgICAgICAgcmV0dXJuIFwiXCI7XG4gICAgfVxuICAgIGlmIChub2RlLmJ2aE5hbWUgPT09IFwiaGlwXCIpIHtcbiAgICAgICAgcmV0dXJuIFwiQ0hBTk5FTFMgNiBYcG9zaXRpb24gWXBvc2l0aW9uIFpwb3NpdGlvbiBYcm90YXRpb24gWXJvdGF0aW9uIFpyb3RhdGlvblwiO1xuICAgIH1cbiAgICByZXR1cm4gXCJDSEFOTkVMUyBcIiArIG5vZGUuY2hhbm5lbHMubGVuZ3RoICsgXCIgXCIgKyBub2RlLmNoYW5uZWxzLmpvaW4oXCIgXCIpO1xufVxuZnVuY3Rpb24gYXBwZW5kTm9kZShqb2ludCwgdGFicykge1xuICAgIGxldCByZXN1bHQgPSBcIlwiO1xuICAgIGNvbnN0IGJvbmVUeXBlID0gKGpvaW50LmJ2aE5hbWUgPT09IFwiaGlwXCIpID8gXCJST09UXCIgOiBcIkpPSU5UXCI7XG4gICAgY29uc3QgY2hhbm5lbHMgPSBjaGFubmVsc1N0cmluZyhqb2ludCk7XG4gICAgY29uc3Qgb2Zmc2V0ID0gKGpvaW50LmJ2aE5hbWUgPT09IFwiaGlwXCIpID8gb2Zmc2V0VG9TdHJpbmcoam9pbnQub2Zmc2V0LCA2KSA6IG9mZnNldFRvU3RyaW5nKGpvaW50Lm9mZnNldCwgNCk7XG4gICAgaWYgKGpvaW50LmJ2aE5hbWUgIT0gXCJlbmRcIikge1xuICAgICAgICByZXN1bHQgKz0gdGFicyArIGJvbmVUeXBlICsgXCIgXCIgKyBqb2ludC5idmhOYW1lICsgXCJcXG5cIiArIHRhYnMgKyBcIntcXG5cIjtcbiAgICB9XG4gICAgZWxzZSB7XG4gICAgICAgIHJlc3VsdCArPSB0YWJzICsgXCJFbmQgU2l0ZVwiICsgXCJcXG5cIiArIHRhYnMgKyBcIntcXG5cIjtcbiAgICB9XG4gICAgcmVzdWx0ICs9IHRhYnMgKyBcIlxcdFwiICsgb2Zmc2V0ICsgXCJcXG5cIjtcbiAgICBpZiAoam9pbnQuYnZoTmFtZSAhPSBcImVuZFwiKSB7XG4gICAgICAgIHJlc3VsdCArPSB0YWJzICsgXCJcXHRcIiArIGNoYW5uZWxzICsgXCJcXG5cIjtcbiAgICB9XG4gICAgaWYgKGpvaW50LmNoaWxkcmVuKSB7XG4gICAgICAgIGpvaW50LmNoaWxkcmVuLmZvckVhY2goKGl0ZW0pID0+IHsgcmVzdWx0ICs9IGFwcGVuZE5vZGUoaXRlbSwgdGFicyArIFwiXFx0XCIpOyB9KTtcbiAgICB9XG4gICAgcmVzdWx0ICs9IHRhYnMgKyBcIn1cXG5cIjtcbiAgICByZXR1cm4gcmVzdWx0O1xufVxuZnVuY3Rpb24gY29udGFpbnNOYW1lcyhub2RlLCBidmhOYW1lcykge1xuICAgIGlmIChidmhOYW1lcy5pbmNsdWRlcyhub2RlLmJ2aE5hbWUpKSB7XG4gICAgICAgIHJldHVybiB0cnVlO1xuICAgIH1cbiAgICBpZiAoIW5vZGUuY2hpbGRyZW4pIHtcbiAgICAgICAgcmV0dXJuIGZhbHNlO1xuICAgIH1cbiAgICByZXR1cm4gISFub2RlLmNoaWxkcmVuLm1hcCgoaXRlbSkgPT4gY29udGFpbnNOYW1lcyhpdGVtLCBidmhOYW1lcykpLmZpbmQoKGl0ZW0pID0+ICEhaXRlbSk7XG59XG5mdW5jdGlvbiBjb2xsZWN0Tm9kZXMobm9kZSwgYnZoTmFtZXMpIHtcbiAgICBjb25zdCByZXN1bHQgPSB7fTtcbiAgICBpZiAoY29udGFpbnNOYW1lcyhub2RlLCBidmhOYW1lcykpIHtcbiAgICAgICAgcmVzdWx0LmJ2aE5hbWUgPSBub2RlLmJ2aE5hbWU7XG4gICAgfVxuICAgIGVsc2Uge1xuICAgICAgICByZXN1bHQuZXhjbHVkZSA9IHRydWU7XG4gICAgICAgIHJldHVybiByZXN1bHQ7XG4gICAgfVxuICAgIGlmIChub2RlLmNoaWxkcmVuICYmICEhbm9kZS5jaGlsZHJlbi5tYXAoKGl0ZW0pID0+IGNvbnRhaW5zTmFtZXMoaXRlbSwgYnZoTmFtZXMpKS5maW5kKChpdGVtKSA9PiAhIWl0ZW0pKSB7XG4gICAgICAgIHJlc3VsdC5jaGlsZHJlbiA9IG5vZGUuY2hpbGRyZW4ubWFwKChpdGVtKSA9PiBjb2xsZWN0Tm9kZXMoaXRlbSwgYnZoTmFtZXMpKS5maWx0ZXIoKGl0ZW0pID0+ICFpdGVtLmV4Y2x1ZGUpO1xuICAgIH1cbiAgICBlbHNlIHtcbiAgICAgICAgcmVzdWx0LmNoaWxkcmVuID0gW107XG4gICAgfVxuICAgIGlmIChyZXN1bHQuY2hpbGRyZW4ubGVuZ3RoID4gMCkge1xuICAgICAgICByZXR1cm4gcmVzdWx0O1xuICAgIH1cbiAgICByZXN1bHQuY2hpbGRyZW4ucHVzaCh7IGJ2aE5hbWU6IFwiZW5kXCIgfSk7XG4gICAgcmV0dXJuIHJlc3VsdDtcbn1cbmZ1bmN0aW9uIHN1YlRyZWUoam9pbnRzKSB7XG4gICAgY29uc3QgbmFtZXMgPSBqb2ludHMubWFwKGl0ZW0gPT4gaXRlbS5qb2ludF9uYW1lKTtcbiAgICBjb25zdCBidmhOYW1lcyA9IG5hbWVzLm1hcChpdGVtID0+IGFsaWFzZXNfMS5hbGlhc2VzW2l0ZW1dIHx8IGl0ZW0pO1xuICAgIHJldHVybiBjb2xsZWN0Tm9kZXMoaGllcmFyY2h5XzEuaGllcmFyY2h5LCBidmhOYW1lcyk7XG59XG5mdW5jdGlvbiB2aXNpdE5vZGUobm9kZSwgdmlzaXRvciwgY2hpbGRyZW5GaXJzdCA9IGZhbHNlKSB7XG4gICAgaWYgKG5vZGUuY2hpbGRyZW4gJiYgY2hpbGRyZW5GaXJzdCkge1xuICAgICAgICBub2RlLmNoaWxkcmVuLnRvUmV2ZXJzZWQoKS5mb3JFYWNoKChpdGVtKSA9PiB2aXNpdE5vZGUoaXRlbSwgdmlzaXRvciwgdHJ1ZSkpO1xuICAgIH1cbiAgICB2aXNpdG9yKG5vZGUpO1xuICAgIGlmIChub2RlLmNoaWxkcmVuICYmICFjaGlsZHJlbkZpcnN0KSB7XG4gICAgICAgIG5vZGUuY2hpbGRyZW4udG9SZXZlcnNlZCgpLmZvckVhY2goKGl0ZW0pID0+IHZpc2l0Tm9kZShpdGVtLCB2aXNpdG9yLCBmYWxzZSkpO1xuICAgIH1cbn1cbmZ1bmN0aW9uIGV4dHJhY3RGcmFtZXNMZW5ndGgoYW5pbUpvaW50cykge1xuICAgIHZhciBfYSwgX2I7XG4gICAgY29uc3Qgam9pbnQgPSBhbmltSm9pbnRzLmZpbmQoKGl0ZW0pID0+IHsgdmFyIF9hLCBfYjsgcmV0dXJuICgoX2EgPSBpdGVtLnBvc2l0aW9uX2tleXMpID09PSBudWxsIHx8IF9hID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfYS5sZW5ndGgpIHx8ICgoX2IgPSBpdGVtLnJvdGF0aW9uX2tleXMpID09PSBudWxsIHx8IF9iID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfYi5sZW5ndGgpOyB9KTtcbiAgICByZXR1cm4gKChfYSA9IGpvaW50ID09PSBudWxsIHx8IGpvaW50ID09PSB2b2lkIDAgPyB2b2lkIDAgOiBqb2ludC5wb3NpdGlvbl9rZXlzKSA9PT0gbnVsbCB8fCBfYSA9PT0gdm9pZCAwID8gdm9pZCAwIDogX2EubGVuZ3RoKSB8fCAoKF9iID0gam9pbnQgPT09IG51bGwgfHwgam9pbnQgPT09IHZvaWQgMCA/IHZvaWQgMCA6IGpvaW50LnJvdGF0aW9uX2tleXMpID09PSBudWxsIHx8IF9iID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfYi5sZW5ndGgpO1xufVxuZnVuY3Rpb24gZXh0cmFjdFRpbWVzKGFuaW1Kb2ludHMpIHtcbiAgICB2YXIgX2E7XG4gICAgY29uc3Qgam9pbnQgPSBhbmltSm9pbnRzLmZpbmQoKGl0ZW0pID0+IHsgdmFyIF9hLCBfYjsgcmV0dXJuICgoX2EgPSBpdGVtLnBvc2l0aW9uX2tleXMpID09PSBudWxsIHx8IF9hID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfYS5sZW5ndGgpIHx8ICgoX2IgPSBpdGVtLnJvdGF0aW9uX2tleXMpID09PSBudWxsIHx8IF9iID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfYi5sZW5ndGgpOyB9KTtcbiAgICBjb25zdCB0aW1lSG9sZGVycyA9ICgoX2EgPSBqb2ludCA9PT0gbnVsbCB8fCBqb2ludCA9PT0gdm9pZCAwID8gdm9pZCAwIDogam9pbnQucG9zaXRpb25fa2V5cykgPT09IG51bGwgfHwgX2EgPT09IHZvaWQgMCA/IHZvaWQgMCA6IF9hLmxlbmd0aCkgPyBqb2ludC5wb3NpdGlvbl9rZXlzIDogam9pbnQgPT09IG51bGwgfHwgam9pbnQgPT09IHZvaWQgMCA/IHZvaWQgMCA6IGpvaW50LnJvdGF0aW9uX2tleXM7XG4gICAgcmV0dXJuICh0aW1lSG9sZGVycyB8fCBbXSkubWFwKChpdGVtKSA9PiBpdGVtLnRpbWUpO1xufVxuZnVuY3Rpb24gZmlsbENoYW5uZWxzKG5vZGUsIGpvaW50KSB7XG4gICAgdmFyIF9hO1xuICAgIGlmIChub2RlLmJ2aE5hbWUgPT0gXCJoaXBcIikge1xuICAgICAgICBub2RlLmNoYW5uZWxzID0gW1wiWHBvc2l0aW9uXCIsIFwiWXBvc2l0aW9uXCIsIFwiWnBvc2l0aW9uXCIsIFwiWHJvdGF0aW9uXCIsIFwiWXJvdGF0aW9uXCIsIFwiWnJvdGF0aW9uXCJdO1xuICAgICAgICByZXR1cm47XG4gICAgfVxuICAgIG5vZGUuY2hhbm5lbHMgPSBbXTtcbiAgICBpZiAoKF9hID0gam9pbnQgPT09IG51bGwgfHwgam9pbnQgPT09IHZvaWQgMCA/IHZvaWQgMCA6IGpvaW50LnBvc2l0aW9uX2tleXMpID09PSBudWxsIHx8IF9hID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfYS5sZW5ndGgpIHtcbiAgICAgICAgbm9kZS5jaGFubmVscy5wdXNoKFwiWHBvc2l0aW9uXCIpO1xuICAgICAgICBub2RlLmNoYW5uZWxzLnB1c2goXCJZcG9zaXRpb25cIik7XG4gICAgICAgIG5vZGUuY2hhbm5lbHMucHVzaChcIlpwb3NpdGlvblwiKTtcbiAgICB9XG4gICAgbm9kZS5jaGFubmVscy5wdXNoKFwiWHJvdGF0aW9uXCIpO1xuICAgIG5vZGUuY2hhbm5lbHMucHVzaChcIllyb3RhdGlvblwiKTtcbiAgICBub2RlLmNoYW5uZWxzLnB1c2goXCJacm90YXRpb25cIik7XG59XG5mdW5jdGlvbiBhbmltUG9zaXRpb25Ub0J2aChwb3NpdGlvbikge1xuICAgIGNvbnN0IG11bHRpcGxpZXIgPSAzOS4zNzk1O1xuICAgIHJldHVybiB7IHg6IHBvc2l0aW9uLnkgKiBtdWx0aXBsaWVyLCB5OiBwb3NpdGlvbi56ICogbXVsdGlwbGllciwgejogcG9zaXRpb24ueCAqIG11bHRpcGxpZXIgfTtcbn1cbmZ1bmN0aW9uIGZpbGxLZXlGcmFtZXMoZGF0YSwgYnZoTm9kZSwgZnBzKSB7XG4gICAgY29uc3QgYW5pbUpvaW50cyA9IGRhdGEuam9pbnRzO1xuICAgIGNvbnN0IGxlbmd0aCA9IGV4dHJhY3RGcmFtZXNMZW5ndGgoYW5pbUpvaW50cyk7XG4gICAgY29uc3QgYnZoVGltZXMgPSAoMCwgdXRpbHNfMS5nZXRVbmlmb3JtVGltZXMpKGRhdGEuZHVyYXRpb24gfHwgMC41LCAxIC8gZnBzKTtcbiAgICB2aXNpdE5vZGUoYnZoTm9kZSwgKG5vZGUpID0+IHtcbiAgICAgICAgY29uc3Qgam9pbnQgPSBhbmltSm9pbnRzLmZpbmQoKGl0ZW0pID0+IGFsaWFzZXNfMS5hbGlhc2VzW2l0ZW0uam9pbnRfbmFtZV0gPT09IG5vZGUuYnZoTmFtZSk7XG4gICAgICAgIG5vZGUub2Zmc2V0ID0geyB4OiAwLCB5OiAwLCB6OiAwIH07XG4gICAgICAgIGlmIChub2RlLmJ2aE5hbWUgIT0gXCJlbmRcIikge1xuICAgICAgICAgICAgZmlsbENoYW5uZWxzKG5vZGUsIGpvaW50KTtcbiAgICAgICAgICAgIG5vZGUuYW5pbUtleXMgPSB7XG4gICAgICAgICAgICAgICAgcG9zaXRpb25zOiAoam9pbnQgPT09IG51bGwgfHwgam9pbnQgPT09IHZvaWQgMCA/IHZvaWQgMCA6IGpvaW50LnBvc2l0aW9uX2tleXMpIHx8IFtdLFxuICAgICAgICAgICAgICAgIHJvdGF0aW9uczogKGpvaW50ID09PSBudWxsIHx8IGpvaW50ID09PSB2b2lkIDAgPyB2b2lkIDAgOiBqb2ludC5yb3RhdGlvbl9rZXlzKSB8fCBbXVxuICAgICAgICAgICAgfTtcbiAgICAgICAgICAgIGNvbnN0IHBvc2l0aW9uVGltZXMgPSAoMCwgdXRpbHNfMS5jbGlwVGltZXNUb0Nsb3Nlc3RCVkhUaW1lKShub2RlLmFuaW1LZXlzLnBvc2l0aW9ucy5tYXAoKGl0ZW0pID0+IGl0ZW0udGltZSksIGJ2aFRpbWVzKTtcbiAgICAgICAgICAgIGNvbnN0IHJvdGF0aW9uVGltZXMgPSAoMCwgdXRpbHNfMS5jbGlwVGltZXNUb0Nsb3Nlc3RCVkhUaW1lKShub2RlLmFuaW1LZXlzLnJvdGF0aW9ucy5tYXAoKGl0ZW0pID0+IGl0ZW0udGltZSksIGJ2aFRpbWVzKTtcbiAgICAgICAgICAgIGNvbnN0IHBvc2l0aW9ucyA9ICgwLCB1dGlsc18xLmxlcnBWYWx1ZXMpKG5vZGUuYW5pbUtleXMucG9zaXRpb25zLm1hcCgoaXRlbSkgPT4gYW5pbVBvc2l0aW9uVG9CdmgoaXRlbSkpLCBwb3NpdGlvblRpbWVzLCBidmhUaW1lcywgeyB4OiAwLCB5OiAwLCB6OiAwIH0sIHV0aWxzXzEubGVycFZlY3Rvcik7XG4gICAgICAgICAgICBjb25zdCByb3RhdGlvbnMgPSAoMCwgdXRpbHNfMS5sZXJwVmFsdWVzKShub2RlLmFuaW1LZXlzLnJvdGF0aW9ucy5tYXAoKGl0ZW0pID0+ICgwLCB1dGlsc18xLnRvUXVhdGVybmlvbikoaXRlbSkpLCByb3RhdGlvblRpbWVzLCBidmhUaW1lcywgbmV3IHF1YXRlcm5pb25fMS5RdWF0ZXJuaW9uKCksIHV0aWxzXzEubGVycFF1YXRlcm5pb24pLm1hcChpdGVtID0+ICgwLCB1dGlsc18xLnF1YXRlcm5pb25Ub0V1bGVycykoaXRlbSkpO1xuICAgICAgICAgICAgbm9kZS5idmhGcmFtZXMgPSBbXTtcbiAgICAgICAgICAgIGJ2aFRpbWVzLmZvckVhY2goKGl0ZW0sIGkpID0+IG5vZGUuYnZoRnJhbWVzLnB1c2goe1xuICAgICAgICAgICAgICAgIHBvc2l0aW9uOiBwb3NpdGlvbnNbaV0sXG4gICAgICAgICAgICAgICAgcm90YXRpb246IHJvdGF0aW9uc1tpXVxuICAgICAgICAgICAgfSkpO1xuICAgICAgICB9XG4gICAgfSk7XG4gICAgYnZoTm9kZS5idmhUaW1lcyA9IGJ2aFRpbWVzO1xufVxuZnVuY3Rpb24gZ2V0VmFsdWUoYnZoTm9kZSwgY2hhbm5lbCwgZnJhbWVOdW0pIHtcbiAgICBjb25zdCBmcmFtZSA9IGJ2aE5vZGUuYnZoRnJhbWVzW2ZyYW1lTnVtXTtcbiAgICBjb25zdCBrZXkgPSBjaGFubmVsLnRvTG93ZXJDYXNlKClbMF07XG4gICAgY29uc3QgZGF0YSA9IChjaGFubmVsLmluY2x1ZGVzKFwicG9zXCIpKSA/IGZyYW1lLnBvc2l0aW9uIDogZnJhbWUucm90YXRpb247XG4gICAgY29uc3QgdmFsdWUgPSBkYXRhW2tleV07XG4gICAgcmV0dXJuIChNYXRoLmFicyh2YWx1ZSkgPiAwLjAwMDAwMDAxKSA/IHZhbHVlIDogMDtcbn1cbmZ1bmN0aW9uIGdldFZhbHVlcyhidmhOb2RlLCBmcmFtZU51bSkge1xuICAgIHJldHVybiBidmhOb2RlLmNoYW5uZWxzLm1hcCgoaXRlbSkgPT4gZ2V0VmFsdWUoYnZoTm9kZSwgaXRlbSwgZnJhbWVOdW0pKTtcbn1cbmZ1bmN0aW9uIGdldEZyYW1lVmFsdWVzKGJ2aE5vZGUsIGZyYW1lTnVtKSB7XG4gICAgY29uc3QgcmVzdWx0ID0gW107XG4gICAgdmlzaXROb2RlKGJ2aE5vZGUsIChub2RlKSA9PiB7XG4gICAgICAgIGlmICghbm9kZS5jaGFubmVscykge1xuICAgICAgICAgICAgcmV0dXJuO1xuICAgICAgICB9XG4gICAgICAgIHJlc3VsdC51bnNoaWZ0KC4uLmdldFZhbHVlcyhub2RlLCBmcmFtZU51bSkpO1xuICAgIH0sIHRydWUpO1xuICAgIHJldHVybiByZXN1bHQ7XG59XG5mdW5jdGlvbiBnZXRGcmFtZVJvdyhidmhOb2RlLCBmcmFtZU51bSkge1xuICAgIGNvbnN0IHZhbHVlcyA9IGdldEZyYW1lVmFsdWVzKGJ2aE5vZGUsIGZyYW1lTnVtKTtcbiAgICByZXR1cm4gdmFsdWVzLm1hcChpdGVtID0+ICgwLCB1dGlsc18xLmZsb2F0VG9TdHJpbmcpKGl0ZW0sIDQpKS5qb2luKFwiIFwiKSArIFwiIFxcblwiO1xufVxuZnVuY3Rpb24gc2VyaWFsaXplQlZIKGJ2aE5vZGUpIHtcbiAgICBsZXQgcmVzdWx0ID0gXCJISUVSQVJDSFlcXG5cIjtcbiAgICByZXN1bHQgKz0gYXBwZW5kTm9kZShidmhOb2RlLCBcIlwiKTtcbiAgICByZXN1bHQgKz0gXCJNT1RJT05cXG5cIjtcbiAgICByZXN1bHQgKz0gXCJGcmFtZXM6IFwiICsgYnZoTm9kZS5idmhUaW1lcy5sZW5ndGggKyBcIlxcblwiO1xuICAgIHJlc3VsdCArPSBcIkZyYW1lIFRpbWU6IFwiICsgKDAsIHV0aWxzXzEuZmxvYXRUb1N0cmluZykoYnZoTm9kZS5idmhUaW1lc1sxXSwgNikgKyBcIlxcblwiO1xuICAgIGZvciAobGV0IGkgPSAwOyBpIDwgYnZoTm9kZS5idmhUaW1lcy5sZW5ndGg7IGkrKykge1xuICAgICAgICByZXN1bHQgKz0gZ2V0RnJhbWVSb3coYnZoTm9kZSwgaSk7XG4gICAgfVxuICAgIHJldHVybiByZXN1bHQ7XG59XG5mdW5jdGlvbiB0b0JWSChkYXRhLCBmcHMgPSAyNCkge1xuICAgIGNvbnN0IGJ2aE5vZGUgPSBzdWJUcmVlKGRhdGEuam9pbnRzKTtcbiAgICBmaWxsS2V5RnJhbWVzKGRhdGEsIGJ2aE5vZGUsIGZwcyk7XG4gICAgcmV0dXJuIGJ2aE5vZGU7XG59XG4iLCJcInVzZSBzdHJpY3RcIjtcbk9iamVjdC5kZWZpbmVQcm9wZXJ0eShleHBvcnRzLCBcIl9fZXNNb2R1bGVcIiwgeyB2YWx1ZTogdHJ1ZSB9KTtcbmV4cG9ydHMuaGllcmFyY2h5ID0gdm9pZCAwO1xuZXhwb3J0cy5oaWVyYXJjaHkgPSB7XG4gICAgXCJidmhOYW1lXCI6IFwiaGlwXCIsXG4gICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgIHtcbiAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1TcGluZTFcIixcbiAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibVNwaW5lMlwiLFxuICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJhYmRvbWVuXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1TcGluZTNcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibVNwaW5lNFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJjaGVzdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJuZWNrXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImhlYWRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZmlndXJlSGFpclwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRXllUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUV5ZUxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VSb290XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlRXllQWx0UmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlRXllQWx0TGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VGb3JlaGVhZExlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlRm9yZWhlYWRSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VFeWVicm93T3V0ZXJMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZUV5ZWJyb3dDZW50ZXJMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZUV5ZWJyb3dJbm5lckxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlRXllYnJvd091dGVyUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlRXllYnJvd0NlbnRlclJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZUV5ZWJyb3dJbm5lclJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZUV5ZUxpZFVwcGVyTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VFeWVMaWRMb3dlckxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlRXllTGlkVXBwZXJSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VFeWVMaWRMb3dlclJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZUVhcjFMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlRWFyMkxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VFYXIxUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VFYXIyUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VOb3NlTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VOb3NlQ2VudGVyXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZU5vc2VSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VDaGVla0xvd2VyTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VDaGVla1VwcGVyTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VDaGVla0xvd2VyUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlQ2hlZWtVcHBlclJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZUphd1wiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZUNoaW5cIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZVRlZXRoTG93ZXJcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VMaXBMb3dlckxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VMaXBMb3dlclJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlTGlwTG93ZXJDZW50ZXJcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VUb25ndWVCYXNlXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlVG9uZ3VlVGlwXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZUphd1NoYXBlclwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VGb3JlaGVhZENlbnRlclwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VOb3NlQmFzZVwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VUZWV0aFVwcGVyXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlTGlwVXBwZXJMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VMaXBVcHBlclJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VMaXBDb3JuZXJMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VMaXBDb3JuZXJSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlTGlwVXBwZXJDZW50ZXJcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VFeWVjb3JuZXJJbm5lckxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlRXllY29ybmVySW5uZXJSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VOb3NlQnJpZGdlXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJsQ29sbGFyXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImxTaGxkclwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJsRm9yZUFybVwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJsSGFuZFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZE1pZGRsZTFMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1IYW5kTWlkZGxlMkxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRNaWRkbGUzTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZEluZGV4MUxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRJbmRleDJMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1IYW5kSW5kZXgzTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZFJpbmcxTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZFJpbmcyTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZFJpbmczTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZFBpbmt5MUxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRQaW5reTJMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1IYW5kUGlua3kzTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZFRodW1iMUxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRUaHVtYjJMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1IYW5kVGh1bWIzTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwickNvbGxhclwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJyU2hsZHJcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwickZvcmVBcm1cIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwickhhbmRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRNaWRkbGUxUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRNaWRkbGUyUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRNaWRkbGUzUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRJbmRleDFSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZEluZGV4MlJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1IYW5kSW5kZXgzUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRSaW5nMVJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1IYW5kUmluZzJSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZFJpbmczUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRQaW5reTFSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZFBpbmt5MlJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1IYW5kUGlua3kzUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRUaHVtYjFSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZFRodW1iMlJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1IYW5kVGh1bWIzUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1XaW5nc1Jvb3RcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibVdpbmcxTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtV2luZzJMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1XaW5nM0xlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibVdpbmc0TGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1XaW5nNEZhbkxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtV2luZzFSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtV2luZzJSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtV2luZzNSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtV2luZzRSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1XaW5nNEZhblJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICBdXG4gICAgICAgIH0sXG4gICAgICAgIHtcbiAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcInJUaGlnaFwiLFxuICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJyU2hpblwiLFxuICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJyRm9vdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRm9vdFJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1Ub2VSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgXVxuICAgICAgICB9LFxuICAgICAgICB7XG4gICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJsVGhpZ2hcIixcbiAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibFNoaW5cIixcbiAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibEZvb3RcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZvb3RMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1Ub2VMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICBdXG4gICAgICAgIH0sXG4gICAgICAgIHtcbiAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1UYWlsMVwiLFxuICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtVGFpbDJcIixcbiAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibVRhaWwzXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1UYWlsNFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtVGFpbDVcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibVRhaWw2XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICBdXG4gICAgICAgIH0sXG4gICAgICAgIHtcbiAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1Hcm9pblwiLFxuICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIF1cbiAgICAgICAgfSxcbiAgICAgICAge1xuICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhpbmRMaW1ic1Jvb3RcIixcbiAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhpbmRMaW1iMUxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhpbmRMaW1iMkxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhpbmRMaW1iM0xlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhpbmRMaW1iNExlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGluZExpbWIxUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhpbmRMaW1iMlJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1IaW5kTGltYjNSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGluZExpbWI0UmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIF1cbiAgICAgICAgfVxuICAgIF1cbn07XG4iLCJcInVzZSBzdHJpY3RcIjtcbk9iamVjdC5kZWZpbmVQcm9wZXJ0eShleHBvcnRzLCBcIl9fZXNNb2R1bGVcIiwgeyB2YWx1ZTogdHJ1ZSB9KTtcbmV4cG9ydHMubWFsZU9mZnNldHMgPSBleHBvcnRzLmZlbWFsZU9mZnNldHMgPSB2b2lkIDA7XG5leHBvcnRzLmZlbWFsZU9mZnNldHMgPSB7XG4gICAgXCJoaXBcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcIm1IaW5kTGltYnNSb290XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAzLjMwNzEsXG4gICAgICAgIFwielwiOiAtNy44NzRcbiAgICB9LFxuICAgIFwibUhpbmRMaW1iMVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC01LjA3ODcsXG4gICAgICAgIFwieVwiOiAtNC45MjEzLFxuICAgICAgICBcInpcIjogLTguMDMxNVxuICAgIH0sXG4gICAgXCJtSGluZExpbWIyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMS44MTEsXG4gICAgICAgIFwieVwiOiAtMTkuMzMwNyxcbiAgICAgICAgXCJ6XCI6IDAuMDc4N1xuICAgIH0sXG4gICAgXCJtSGluZExpbWIzUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMC4xMTgxLFxuICAgICAgICBcInlcIjogLTE4LjQyNTIsXG4gICAgICAgIFwielwiOiAtMS4xODExXG4gICAgfSxcbiAgICBcIm1IaW5kTGltYjRSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTIuNDAxNixcbiAgICAgICAgXCJ6XCI6IDQuNDA5NFxuICAgIH0sXG4gICAgXCJlbmRfbUhpbmRMaW1iNFJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjMxNSxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiA0LjEzMzlcbiAgICB9LFxuICAgIFwibUhpbmRMaW1iMUxlZnRcIjoge1xuICAgICAgICBcInhcIjogNS4wNzg3LFxuICAgICAgICBcInlcIjogLTQuOTIxMyxcbiAgICAgICAgXCJ6XCI6IC04LjAzMTVcbiAgICB9LFxuICAgIFwibUhpbmRMaW1iMkxlZnRcIjoge1xuICAgICAgICBcInhcIjogLTEuODExLFxuICAgICAgICBcInlcIjogLTE5LjMzMDcsXG4gICAgICAgIFwielwiOiAwLjA3ODdcbiAgICB9LFxuICAgIFwibUhpbmRMaW1iM0xlZnRcIjoge1xuICAgICAgICBcInhcIjogLTAuMTE4MSxcbiAgICAgICAgXCJ5XCI6IC0xOC40MjUyLFxuICAgICAgICBcInpcIjogLTEuMTgxMVxuICAgIH0sXG4gICAgXCJtSGluZExpbWI0TGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTIuNDAxNixcbiAgICAgICAgXCJ6XCI6IDQuNDA5NFxuICAgIH0sXG4gICAgXCJlbmRfbUhpbmRMaW1iNExlZnRcIjoge1xuICAgICAgICBcInhcIjogMC4zMTUsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogNC4xMzM5XG4gICAgfSxcbiAgICBcIm1Hcm9pblwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTMuODE4OSxcbiAgICAgICAgXCJ6XCI6IDIuNTE5N1xuICAgIH0sXG4gICAgXCJlbmRfbUdyb2luXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMi41OTg0LFxuICAgICAgICBcInpcIjogMC4xNTc1XG4gICAgfSxcbiAgICBcIm1UYWlsMVwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMS44NTA0LFxuICAgICAgICBcInpcIjogLTQuNTY2OVxuICAgIH0sXG4gICAgXCJtVGFpbDJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtNy43NTU5XG4gICAgfSxcbiAgICBcIm1UYWlsM1wiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC02LjYxNDJcbiAgICB9LFxuICAgIFwibVRhaWw0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTUuNTkwNVxuICAgIH0sXG4gICAgXCJtVGFpbDVcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtNC40MDk0XG4gICAgfSxcbiAgICBcIm1UYWlsNlwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC0zLjcwMDhcbiAgICB9LFxuICAgIFwiZW5kX21UYWlsNlwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC0zLjUwMzlcbiAgICB9LFxuICAgIFwibFRoaWdoXCI6IHtcbiAgICAgICAgXCJ4XCI6IDQuOTkwNyxcbiAgICAgICAgXCJ5XCI6IC0xLjYxNDEsXG4gICAgICAgIFwielwiOiAxLjMyOVxuICAgIH0sXG4gICAgXCJsU2hpblwiOiB7XG4gICAgICAgIFwieFwiOiAtMS43OTQsXG4gICAgICAgIFwieVwiOiAtMTkuMzMyOCxcbiAgICAgICAgXCJ6XCI6IC0wLjAzNDlcbiAgICB9LFxuICAgIFwibEZvb3RcIjoge1xuICAgICAgICBcInhcIjogMC4wNTQzLFxuICAgICAgICBcInlcIjogLTE4LjQ0MjgsXG4gICAgICAgIFwielwiOiAtMS4xMzczXG4gICAgfSxcbiAgICBcIm1Gb290TGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTIuMzg2NixcbiAgICAgICAgXCJ6XCI6IDQuNDA3N1xuICAgIH0sXG4gICAgXCJtVG9lTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDQuMjkxM1xuICAgIH0sXG4gICAgXCJlbmRfbVRvZUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjc4NzRcbiAgICB9LFxuICAgIFwiclRoaWdoXCI6IHtcbiAgICAgICAgXCJ4XCI6IC01LjA3MTEsXG4gICAgICAgIFwieVwiOiAtMS42MTc2LFxuICAgICAgICBcInpcIjogMS4zMjM2XG4gICAgfSxcbiAgICBcInJTaGluXCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuOTE0OCxcbiAgICAgICAgXCJ5XCI6IC0xOS4zMjc2LFxuICAgICAgICBcInpcIjogLTAuMDMwN1xuICAgIH0sXG4gICAgXCJyRm9vdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTE4LjQ0NDYsXG4gICAgICAgIFwielwiOiAtMS4xMzY2XG4gICAgfSxcbiAgICBcIm1Gb290UmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0yLjM4NzMsXG4gICAgICAgIFwielwiOiA0LjQwNzdcbiAgICB9LFxuICAgIFwibVRvZVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogNC4yOTEzXG4gICAgfSxcbiAgICBcImVuZF9tVG9lUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjc4NzRcbiAgICB9LFxuICAgIFwibVNwaW5lMVwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMy4zMSxcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwibVNwaW5lMlwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTMuMzEsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcImFiZG9tZW5cIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDMuMzEsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcIm1TcGluZTNcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDguMDY2LFxuICAgICAgICBcInpcIjogLTAuNjA1XG4gICAgfSxcbiAgICBcIm1TcGluZTRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC04LjA2NixcbiAgICAgICAgXCJ6XCI6IDAuNjA1XG4gICAgfSxcbiAgICBcImNoZXN0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiA4LjA2NixcbiAgICAgICAgXCJ6XCI6IC0wLjYwNVxuICAgIH0sXG4gICAgXCJtV2luZ3NSb290XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTAuNTQ1N1xuICAgIH0sXG4gICAgXCJtV2luZzFSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMy43NDAyLFxuICAgICAgICBcInlcIjogNy4xMjYsXG4gICAgICAgIFwielwiOiAtMy44OTc2XG4gICAgfSxcbiAgICBcIm1XaW5nMlJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC02LjY1MzUsXG4gICAgICAgIFwieVwiOiAyLjYzNzgsXG4gICAgICAgIFwielwiOiAtNi42MTQyXG4gICAgfSxcbiAgICBcIm1XaW5nM1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC03LjIwNDcsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTcuMTI2XG4gICAgfSxcbiAgICBcIm1XaW5nNEZhblJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC02LjgxMSxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtNi43MzIzXG4gICAgfSxcbiAgICBcImVuZF9tV2luZzRGYW5SaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMi40NDA5LFxuICAgICAgICBcInlcIjogLTYuMjU5OCxcbiAgICAgICAgXCJ6XCI6IC0yLjY3NzJcbiAgICB9LFxuICAgIFwibVdpbmc0UmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTYuODExLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC02LjczMjNcbiAgICB9LFxuICAgIFwiZW5kX21XaW5nNFJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC01LjE5NjgsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTUuNzQ4XG4gICAgfSxcbiAgICBcIm1XaW5nMUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMy43NDAyLFxuICAgICAgICBcInlcIjogNy4xMjYsXG4gICAgICAgIFwielwiOiAtMy44OTc2XG4gICAgfSxcbiAgICBcIm1XaW5nMkxlZnRcIjoge1xuICAgICAgICBcInhcIjogNi42NTM1LFxuICAgICAgICBcInlcIjogMi42Mzc4LFxuICAgICAgICBcInpcIjogLTYuNjE0MlxuICAgIH0sXG4gICAgXCJtV2luZzNMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDcuMjA0NyxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtNy4xMjZcbiAgICB9LFxuICAgIFwibVdpbmc0RmFuTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiA2LjgxMSxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtNi43MzIzXG4gICAgfSxcbiAgICBcImVuZF9tV2luZzRGYW5MZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDIuNDQwOSxcbiAgICAgICAgXCJ5XCI6IC02LjI1OTgsXG4gICAgICAgIFwielwiOiAtMi42NzcyXG4gICAgfSxcbiAgICBcIm1XaW5nNExlZnRcIjoge1xuICAgICAgICBcInhcIjogNi44MTEsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTYuNzMyM1xuICAgIH0sXG4gICAgXCJlbmRfbVdpbmc0TGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiA1LjE5NjgsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTUuNzQ4XG4gICAgfSxcbiAgICBcInJDb2xsYXJcIjoge1xuICAgICAgICBcInhcIjogLTIuODM0NixcbiAgICAgICAgXCJ5XCI6IDYuNTExNixcbiAgICAgICAgXCJ6XCI6IC0wLjgxNTdcbiAgICB9LFxuICAgIFwiclNobGRyXCI6IHtcbiAgICAgICAgXCJ4XCI6IC0zLjEyNjcsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJyRm9yZUFybVwiOiB7XG4gICAgICAgIFwieFwiOiAtMTAuOTM1NCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcInJIYW5kXCI6IHtcbiAgICAgICAgXCJ4XCI6IC05LjUyMzYsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJtSGFuZFRodW1iMVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjAyMzYsXG4gICAgICAgIFwieVwiOiAwLjE1NzUsXG4gICAgICAgIFwielwiOiAxLjIyMDVcbiAgICB9LFxuICAgIFwibUhhbmRUaHVtYjJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS4yNTk4LFxuICAgICAgICBcInlcIjogLTAuMDM5NCxcbiAgICAgICAgXCJ6XCI6IDEuMTAyNFxuICAgIH0sXG4gICAgXCJtSGFuZFRodW1iM1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjIyMDUsXG4gICAgICAgIFwieVwiOiAtMC4wMzk0LFxuICAgICAgICBcInpcIjogMC45MDU1XG4gICAgfSxcbiAgICBcImVuZF9tSGFuZFRodW1iM1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjk4NDMsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC41OTA2XG4gICAgfSxcbiAgICBcIm1IYW5kUGlua3kxUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTMuNzQwMSxcbiAgICAgICAgXCJ5XCI6IDAuMTE4MSxcbiAgICAgICAgXCJ6XCI6IC0xLjIyMDVcbiAgICB9LFxuICAgIFwibUhhbmRQaW5reTJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC45ODQyLFxuICAgICAgICBcInlcIjogLTAuMjM2MixcbiAgICAgICAgXCJ6XCI6IC0wLjk0NDlcbiAgICB9LFxuICAgIFwibUhhbmRQaW5reTNSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC43MDg3LFxuICAgICAgICBcInlcIjogLTAuMTU3NSxcbiAgICAgICAgXCJ6XCI6IC0wLjU5MDZcbiAgICB9LFxuICAgIFwiZW5kX21IYW5kUGlua3kzUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuNjI5OSxcbiAgICAgICAgXCJ5XCI6IC0wLjE1NzUsXG4gICAgICAgIFwielwiOiAtMC41MTE4XG4gICAgfSxcbiAgICBcIm1IYW5kUmluZzFSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMy44OTc2LFxuICAgICAgICBcInlcIjogMC4zNTQzLFxuICAgICAgICBcInpcIjogLTAuMzkzN1xuICAgIH0sXG4gICAgXCJtSGFuZFJpbmcyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuNDk2MSxcbiAgICAgICAgXCJ5XCI6IC0wLjMxNSxcbiAgICAgICAgXCJ6XCI6IC0wLjUxMThcbiAgICB9LFxuICAgIFwibUhhbmRSaW5nM1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjU3NDgsXG4gICAgICAgIFwieVwiOiAtMC4zNTQzLFxuICAgICAgICBcInpcIjogLTAuNTExOFxuICAgIH0sXG4gICAgXCJlbmRfbUhhbmRSaW5nM1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjEwMjQsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogLTAuMzkzN1xuICAgIH0sXG4gICAgXCJtSGFuZEluZGV4MVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0zLjgxODksXG4gICAgICAgIFwieVwiOiAwLjU5MDUsXG4gICAgICAgIFwielwiOiAxLjQ5NjFcbiAgICB9LFxuICAgIFwibUhhbmRJbmRleDJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS40MTczLFxuICAgICAgICBcInlcIjogLTAuMjM2MixcbiAgICAgICAgXCJ6XCI6IDAuNjY5M1xuICAgIH0sXG4gICAgXCJtSGFuZEluZGV4M1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjI1OTgsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogMC41NTEyXG4gICAgfSxcbiAgICBcImVuZF9tSGFuZEluZGV4M1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjk4NDMsXG4gICAgICAgIFwieVwiOiAtMC4xNTc1LFxuICAgICAgICBcInpcIjogMC40MzMxXG4gICAgfSxcbiAgICBcIm1IYW5kTWlkZGxlMVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0zLjk3NjQsXG4gICAgICAgIFwieVwiOiAwLjU5MDUsXG4gICAgICAgIFwielwiOiAwLjUxMThcbiAgICB9LFxuICAgIFwibUhhbmRNaWRkbGUyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuNTc0OCxcbiAgICAgICAgXCJ5XCI6IC0wLjIzNjIsXG4gICAgICAgIFwielwiOiAtMC4wMzk0XG4gICAgfSxcbiAgICBcIm1IYW5kTWlkZGxlM1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjkyOTEsXG4gICAgICAgIFwieVwiOiAtMC4zMTUsXG4gICAgICAgIFwielwiOiAtMC4wMzk0XG4gICAgfSxcbiAgICBcImVuZF9tSGFuZE1pZGRsZTNSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS4yOTkyLFxuICAgICAgICBcInlcIjogLTAuMjM2MixcbiAgICAgICAgXCJ6XCI6IC0wLjA3ODdcbiAgICB9LFxuICAgIFwibENvbGxhclwiOiB7XG4gICAgICAgIFwieFwiOiAyLjgyMixcbiAgICAgICAgXCJ5XCI6IDYuNTExNixcbiAgICAgICAgXCJ6XCI6IC0wLjgxNTdcbiAgICB9LFxuICAgIFwibFNobGRyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDMuMTEwMixcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcImxGb3JlQXJtXCI6IHtcbiAgICAgICAgXCJ4XCI6IDEwLjkzNTQsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJsSGFuZFwiOiB7XG4gICAgICAgIFwieFwiOiA5LjUxNjUsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJtSGFuZFRodW1iMUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMS4wMjM2LFxuICAgICAgICBcInlcIjogMC4xNTc1LFxuICAgICAgICBcInpcIjogMS4yMjA1XG4gICAgfSxcbiAgICBcIm1IYW5kVGh1bWIyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjI1OTgsXG4gICAgICAgIFwieVwiOiAtMC4wMzk0LFxuICAgICAgICBcInpcIjogMS4xMDI0XG4gICAgfSxcbiAgICBcIm1IYW5kVGh1bWIzTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjIyMDUsXG4gICAgICAgIFwieVwiOiAtMC4wMzk0LFxuICAgICAgICBcInpcIjogMC45MDU1XG4gICAgfSxcbiAgICBcImVuZF9tSGFuZFRodW1iM0xlZnRcIjoge1xuICAgICAgICBcInhcIjogMC45ODQzLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuNTkwNlxuICAgIH0sXG4gICAgXCJtSGFuZFBpbmt5MUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMy43NDAyLFxuICAgICAgICBcInlcIjogMC4xMTgxLFxuICAgICAgICBcInpcIjogLTEuMjIwNVxuICAgIH0sXG4gICAgXCJtSGFuZFBpbmt5MkxlZnRcIjoge1xuICAgICAgICBcInhcIjogMC45ODQyLFxuICAgICAgICBcInlcIjogLTAuMjM2MixcbiAgICAgICAgXCJ6XCI6IC0wLjk0NDlcbiAgICB9LFxuICAgIFwibUhhbmRQaW5reTNMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuNzA4NyxcbiAgICAgICAgXCJ5XCI6IC0wLjE1NzUsXG4gICAgICAgIFwielwiOiAtMC41OTA2XG4gICAgfSxcbiAgICBcImVuZF9tSGFuZFBpbmt5M0xlZnRcIjoge1xuICAgICAgICBcInhcIjogMC42Mjk5LFxuICAgICAgICBcInlcIjogLTAuMTU3NSxcbiAgICAgICAgXCJ6XCI6IC0wLjUxMThcbiAgICB9LFxuICAgIFwibUhhbmRSaW5nMUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMy44OTc2LFxuICAgICAgICBcInlcIjogMC4zNTQzLFxuICAgICAgICBcInpcIjogLTAuMzkzN1xuICAgIH0sXG4gICAgXCJtSGFuZFJpbmcyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjQ5NjEsXG4gICAgICAgIFwieVwiOiAtMC4zMTUsXG4gICAgICAgIFwielwiOiAtMC41MTE4XG4gICAgfSxcbiAgICBcIm1IYW5kUmluZzNMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuNTc0OCxcbiAgICAgICAgXCJ5XCI6IC0wLjM1NDMsXG4gICAgICAgIFwielwiOiAtMC41MTE4XG4gICAgfSxcbiAgICBcImVuZF9tSGFuZFJpbmczTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjEwMjQsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogLTAuMzkzN1xuICAgIH0sXG4gICAgXCJtSGFuZEluZGV4MUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMy44MTg5LFxuICAgICAgICBcInlcIjogMC41OTA1LFxuICAgICAgICBcInpcIjogMS40OTYxXG4gICAgfSxcbiAgICBcIm1IYW5kSW5kZXgyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjQxNzMsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogMC42NjkzXG4gICAgfSxcbiAgICBcIm1IYW5kSW5kZXgzTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjI1OTgsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogMC41NTEyXG4gICAgfSxcbiAgICBcImVuZF9tSGFuZEluZGV4M0xlZnRcIjoge1xuICAgICAgICBcInhcIjogMC45ODQzLFxuICAgICAgICBcInlcIjogLTAuMTU3NSxcbiAgICAgICAgXCJ6XCI6IDAuNDMzMVxuICAgIH0sXG4gICAgXCJtSGFuZE1pZGRsZTFMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDMuOTc2NCxcbiAgICAgICAgXCJ5XCI6IDAuNTkwNSxcbiAgICAgICAgXCJ6XCI6IDAuNTExOFxuICAgIH0sXG4gICAgXCJtSGFuZE1pZGRsZTJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuNTc0OCxcbiAgICAgICAgXCJ5XCI6IC0wLjIzNjIsXG4gICAgICAgIFwielwiOiAtMC4wMzk0XG4gICAgfSxcbiAgICBcIm1IYW5kTWlkZGxlM0xlZnRcIjoge1xuICAgICAgICBcInhcIjogMS45MjkxLFxuICAgICAgICBcInlcIjogLTAuMzE1LFxuICAgICAgICBcInpcIjogLTAuMDM5NFxuICAgIH0sXG4gICAgXCJlbmRfbUhhbmRNaWRkbGUzTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjI5OTIsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogLTAuMDc4N1xuICAgIH0sXG4gICAgXCJuZWNrXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiA5Ljg4NjEsXG4gICAgICAgIFwielwiOiAtMC4zNzA1XG4gICAgfSxcbiAgICBcImhlYWRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDIuOTc3NixcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwibUZhY2VSb290XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAxLjY5MzQsXG4gICAgICAgIFwielwiOiAwLjkxMDRcbiAgICB9LFxuICAgIFwibUZhY2VOb3NlQnJpZGdlXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLjcxMTgsXG4gICAgICAgIFwielwiOiAzLjMxNFxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VOb3NlQnJpZGdlXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLjMxNSxcbiAgICAgICAgXCJ6XCI6IDAuNTkwNlxuICAgIH0sXG4gICAgXCJtRmFjZUV5ZWNvcm5lcklubmVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuNzEzNSxcbiAgICAgICAgXCJ5XCI6IDEuMTEyMixcbiAgICAgICAgXCJ6XCI6IDIuNzMxM1xuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVjb3JuZXJJbm5lclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC42Mjk5XG4gICAgfSxcbiAgICBcIm1GYWNlRXllY29ybmVySW5uZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuNzEzNSxcbiAgICAgICAgXCJ5XCI6IDEuMTEyMixcbiAgICAgICAgXCJ6XCI6IDIuNzMxM1xuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVjb3JuZXJJbm5lckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjYyOTlcbiAgICB9LFxuICAgIFwibUZhY2VUZWV0aFVwcGVyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMS4xNDc4LFxuICAgICAgICBcInpcIjogMC43MjgzXG4gICAgfSxcbiAgICBcIm1GYWNlTGlwVXBwZXJDZW50ZXJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0wLjA5NjQsXG4gICAgICAgIFwielwiOiAxLjYzODhcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlTGlwVXBwZXJDZW50ZXJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAuMDc4NyxcbiAgICAgICAgXCJ6XCI6IDEuNjkyOVxuICAgIH0sXG4gICAgXCJtRmFjZUxpcENvcm5lclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuNjkxOSxcbiAgICAgICAgXCJ5XCI6IC0wLjM2NDIsXG4gICAgICAgIFwielwiOiAxLjAxOTdcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlTGlwQ29ybmVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTIuMDA3OSxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAxLjc3MTdcbiAgICB9LFxuICAgIFwibUZhY2VMaXBDb3JuZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjY5MTksXG4gICAgICAgIFwieVwiOiAtMC4zNjQyLFxuICAgICAgICBcInpcIjogMS4wMTk3XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUxpcENvcm5lckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMi4wMDc5LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDEuNzcxN1xuICAgIH0sXG4gICAgXCJtRmFjZUxpcFVwcGVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0wLjE0NzgsXG4gICAgICAgIFwielwiOiAxLjYzODhcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlTGlwVXBwZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC41OTA2LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDEuNjE0MlxuICAgIH0sXG4gICAgXCJtRmFjZUxpcFVwcGVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTAuMTQ3OCxcbiAgICAgICAgXCJ6XCI6IDEuNjM4OFxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VMaXBVcHBlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMC41OTA2LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDEuNjE0MlxuICAgIH0sXG4gICAgXCJtRmFjZU5vc2VCYXNlXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMC43MDQ3LFxuICAgICAgICBcInpcIjogMy40MjMyXG4gICAgfSxcbiAgICBcImVuZF9tRmFjZU5vc2VCYXNlXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC41NTEyXG4gICAgfSxcbiAgICBcIm1GYWNlRm9yZWhlYWRDZW50ZXJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDIuMzAyNyxcbiAgICAgICAgXCJ6XCI6IDIuNTIzN1xuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VGb3JlaGVhZENlbnRlclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDEuNDE3M1xuICAgIH0sXG4gICAgXCJtRmFjZUphd1NoYXBlclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlSmF3U2hhcGVyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTAuNjY5M1xuICAgIH0sXG4gICAgXCJtRmFjZUphd1wiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTAuNTMzOSxcbiAgICAgICAgXCJ6XCI6IC0wLjAzNjRcbiAgICB9LFxuICAgIFwibUZhY2VUZWV0aExvd2VyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMS40ODY5LFxuICAgICAgICBcInpcIjogMC43NjQ4XG4gICAgfSxcbiAgICBcIm1GYWNlVG9uZ3VlQmFzZVwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMC4xODIxLFxuICAgICAgICBcInpcIjogMS40MjAzXG4gICAgfSxcbiAgICBcIm1GYWNlVG9uZ3VlVGlwXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLjI1NDksXG4gICAgICAgIFwielwiOiAwLjgwMTJcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlVG9uZ3VlVGlwXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC4zOTM3XG4gICAgfSxcbiAgICBcIm1GYWNlTGlwTG93ZXJDZW50ZXJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAxLjYzODhcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlTGlwTG93ZXJDZW50ZXJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAuMDc4NyxcbiAgICAgICAgXCJ6XCI6IDEuNTc0OFxuICAgIH0sXG4gICAgXCJtRmFjZUxpcExvd2VyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAxLjYzODhcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlTGlwTG93ZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC42NjkzLFxuICAgICAgICBcInlcIjogMC4xOTY5LFxuICAgICAgICBcInpcIjogMS4zMzg2XG4gICAgfSxcbiAgICBcIm1GYWNlTGlwTG93ZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMS42Mzg4XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUxpcExvd2VyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjY2OTMsXG4gICAgICAgIFwieVwiOiAwLjE5NjksXG4gICAgICAgIFwielwiOiAxLjMzODZcbiAgICB9LFxuICAgIFwibUZhY2VDaGluXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMi4wMjEyLFxuICAgICAgICBcInpcIjogMi41MzFcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlQ2hpblwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTAuNzA4NyxcbiAgICAgICAgXCJ6XCI6IDAuODI2OFxuICAgIH0sXG4gICAgXCJtRmFjZUNoZWVrVXBwZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS4yNjYzLFxuICAgICAgICBcInlcIjogLTAuMjMxMyxcbiAgICAgICAgXCJ6XCI6IDIuNTQ5MlxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VDaGVla1VwcGVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuNTkwNixcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjg2NjFcbiAgICB9LFxuICAgIFwibUZhY2VDaGVla0xvd2VyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuMjY2MyxcbiAgICAgICAgXCJ5XCI6IC0xLjE5NDEsXG4gICAgICAgIFwielwiOiAxLjgyMDlcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlQ2hlZWtMb3dlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjE4MTEsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC41MTE4XG4gICAgfSxcbiAgICBcIm1GYWNlQ2hlZWtVcHBlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMS4yNjYzLFxuICAgICAgICBcInlcIjogLTAuMjMxMyxcbiAgICAgICAgXCJ6XCI6IDIuNTQ5MlxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VDaGVla1VwcGVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjU5MDYsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC44NjYxXG4gICAgfSxcbiAgICBcIm1GYWNlQ2hlZWtMb3dlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMS4yNjYzLFxuICAgICAgICBcInlcIjogLTEuMTk0MSxcbiAgICAgICAgXCJ6XCI6IDEuODIwOVxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VDaGVla0xvd2VyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjE4MTEsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC41MTE4XG4gICAgfSxcbiAgICBcIm1GYWNlTm9zZVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjU1ODcsXG4gICAgICAgIFwieVwiOiAtMC4xOTU3LFxuICAgICAgICBcInpcIjogMy4xMzE5XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZU5vc2VSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC4xNTc1LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuNTkwNlxuICAgIH0sXG4gICAgXCJtRmFjZU5vc2VDZW50ZXJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0wLjA1MzQsXG4gICAgICAgIFwielwiOiAzLjcxNDZcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlTm9zZUNlbnRlclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuOTg0M1xuICAgIH0sXG4gICAgXCJtRmFjZU5vc2VMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuNTU4NyxcbiAgICAgICAgXCJ5XCI6IC0wLjE5NTcsXG4gICAgICAgIFwielwiOiAzLjEzMTlcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlTm9zZUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMC4xNTc1LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuNTkwNlxuICAgIH0sXG4gICAgXCJtRmFjZUVhcjFSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMi45Nzk1LFxuICAgICAgICBcInlcIjogMC4wNzEyLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJtRmFjZUVhcjJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC43MDg3LFxuICAgICAgICBcInlcIjogMC45ODQyLFxuICAgICAgICBcInpcIjogLTAuNzQ4XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUVhcjJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMS4yOTkyLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJtRmFjZUVhcjFMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDIuOTc5NSxcbiAgICAgICAgXCJ5XCI6IDAuMDcxMixcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwibUZhY2VFYXIyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjcwODcsXG4gICAgICAgIFwieVwiOiAwLjk4NDIsXG4gICAgICAgIFwielwiOiAtMC43NDhcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRWFyMkxlZnRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDEuMjk5MixcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwibUZhY2VFeWVMaWRMb3dlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjQxNTYsXG4gICAgICAgIFwieVwiOiAxLjE4MzQsXG4gICAgICAgIFwielwiOiAyLjY1ODVcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRXllTGlkTG93ZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTAuMjc1NixcbiAgICAgICAgXCJ6XCI6IDAuOTQ0OVxuICAgIH0sXG4gICAgXCJtRmFjZUV5ZUxpZFVwcGVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuNDE1NixcbiAgICAgICAgXCJ5XCI6IDEuMTg4NyxcbiAgICAgICAgXCJ6XCI6IDIuNjU4NVxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVMaWRVcHBlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLjE5NjksXG4gICAgICAgIFwielwiOiAxLjA2M1xuICAgIH0sXG4gICAgXCJtRmFjZUV5ZUxpZExvd2VyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjQxNTYsXG4gICAgICAgIFwieVwiOiAxLjE4MzQsXG4gICAgICAgIFwielwiOiAyLjY1ODVcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRXllTGlkTG93ZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMC4yNzU2LFxuICAgICAgICBcInpcIjogMC45NDQ5XG4gICAgfSxcbiAgICBcIm1GYWNlRXllTGlkVXBwZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuNDE1NixcbiAgICAgICAgXCJ5XCI6IDEuMTg4NyxcbiAgICAgICAgXCJ6XCI6IDIuNjU4NVxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVMaWRVcHBlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAuMTk2OSxcbiAgICAgICAgXCJ6XCI6IDEuMDYzXG4gICAgfSxcbiAgICBcIm1GYWNlRXllYnJvd0lubmVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuOTA2LFxuICAgICAgICBcInlcIjogMS44NTc4LFxuICAgICAgICBcInpcIjogMi43MzEzXG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUV5ZWJyb3dJbm5lclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMS4wMjM2XG4gICAgfSxcbiAgICBcIm1GYWNlRXllYnJvd0NlbnRlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjY3MTEsXG4gICAgICAgIFwieVwiOiAyLjE0NjcsXG4gICAgICAgIFwielwiOiAyLjU0OTJcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRXllYnJvd0NlbnRlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMS4wNjNcbiAgICB9LFxuICAgIFwibUZhY2VFeWVicm93T3V0ZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS45MTY4LFxuICAgICAgICBcInlcIjogMS43MjksXG4gICAgICAgIFwielwiOiAyLjMzMDdcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRXllYnJvd091dGVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuNTExOCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjkwNTVcbiAgICB9LFxuICAgIFwibUZhY2VFeWVicm93SW5uZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuOTA2LFxuICAgICAgICBcInlcIjogMS44NTc4LFxuICAgICAgICBcInpcIjogMi43MzEzXG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUV5ZWJyb3dJbm5lckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAxLjAyMzZcbiAgICB9LFxuICAgIFwibUZhY2VFeWVicm93Q2VudGVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjY3MTEsXG4gICAgICAgIFwieVwiOiAyLjE0NjcsXG4gICAgICAgIFwielwiOiAyLjU0OTJcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRXllYnJvd0NlbnRlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAxLjA2M1xuICAgIH0sXG4gICAgXCJtRmFjZUV5ZWJyb3dPdXRlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMS45MTY4LFxuICAgICAgICBcInlcIjogMS43MjksXG4gICAgICAgIFwielwiOiAyLjMzMDdcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRXllYnJvd091dGVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjUxMTgsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC45MDU1XG4gICAgfSxcbiAgICBcIm1GYWNlRm9yZWhlYWRSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS4zMDM1LFxuICAgICAgICBcInlcIjogMy4wNjA4LFxuICAgICAgICBcInpcIjogMi4zMzA3XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUZvcmVoZWFkUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuMTU3NSxcbiAgICAgICAgXCJ5XCI6IDAuNzA4NyxcbiAgICAgICAgXCJ6XCI6IDAuOTQ0OVxuICAgIH0sXG4gICAgXCJtRmFjZUZvcmVoZWFkTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjMwMzUsXG4gICAgICAgIFwieVwiOiAzLjA2MDgsXG4gICAgICAgIFwielwiOiAyLjMzMDdcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRm9yZWhlYWRMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuMTU3NSxcbiAgICAgICAgXCJ5XCI6IDAuNzA4NyxcbiAgICAgICAgXCJ6XCI6IDAuOTQ0OVxuICAgIH0sXG4gICAgXCJtRmFjZUV5ZUFsdExlZnRcIjoge1xuICAgICAgICBcInhcIjogMS40MTU2LFxuICAgICAgICBcInlcIjogMS4xODM2LFxuICAgICAgICBcInpcIjogMi42NzUyXG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUV5ZUFsdExlZnRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjk4NDNcbiAgICB9LFxuICAgIFwibUZhY2VFeWVBbHRSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS40MTU2LFxuICAgICAgICBcInlcIjogMS4xODM2LFxuICAgICAgICBcInpcIjogMi42NzU0XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUV5ZUFsdFJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC45ODQzXG4gICAgfSxcbiAgICBcIm1FeWVMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuNDIwMyxcbiAgICAgICAgXCJ5XCI6IDIuODc3LFxuICAgICAgICBcInpcIjogMy41ODU3XG4gICAgfSxcbiAgICBcImVuZF9tRXllTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuOTg0M1xuICAgIH0sXG4gICAgXCJtRXllUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuNDIwMyxcbiAgICAgICAgXCJ5XCI6IDIuODc3LFxuICAgICAgICBcInpcIjogMy41ODU5XG4gICAgfSxcbiAgICBcImVuZF9tRXllUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjk4NDNcbiAgICB9LFxuICAgIFwiZmlndXJlSGFpclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMi42MDM4LFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJlbmRfZmlndXJlSGFpclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMS4yOTkyLFxuICAgICAgICBcInpcIjogMFxuICAgIH1cbn07XG5leHBvcnRzLm1hbGVPZmZzZXRzID0ge1xuICAgIFwiaGlwXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJtSGluZExpbWJzUm9vdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMy4zMDcxLFxuICAgICAgICBcInpcIjogLTcuODc0XG4gICAgfSxcbiAgICBcIm1IaW5kTGltYjFSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtNS4wNzg3LFxuICAgICAgICBcInlcIjogLTQuOTIxMyxcbiAgICAgICAgXCJ6XCI6IC04LjAzMTVcbiAgICB9LFxuICAgIFwibUhpbmRMaW1iMlJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuOTAxNixcbiAgICAgICAgXCJ5XCI6IC0xOS4zMzA3LFxuICAgICAgICBcInpcIjogMC4wODI3XG4gICAgfSxcbiAgICBcIm1IaW5kTGltYjNSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjEyNCxcbiAgICAgICAgXCJ5XCI6IC0yMC4yNjc3LFxuICAgICAgICBcInpcIjogLTEuMjQwMlxuICAgIH0sXG4gICAgXCJtSGluZExpbWI0UmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0yLjQwMTYsXG4gICAgICAgIFwielwiOiA0LjQwOTRcbiAgICB9LFxuICAgIFwiZW5kX21IaW5kTGltYjRSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC4zMTUsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogNC4xMzM5XG4gICAgfSxcbiAgICBcIm1IaW5kTGltYjFMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDUuMDc4NyxcbiAgICAgICAgXCJ5XCI6IC00LjkyMTMsXG4gICAgICAgIFwielwiOiAtOC4wMzE1XG4gICAgfSxcbiAgICBcIm1IaW5kTGltYjJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjkwMTYsXG4gICAgICAgIFwieVwiOiAtMTkuMzMwNyxcbiAgICAgICAgXCJ6XCI6IDAuMDgyN1xuICAgIH0sXG4gICAgXCJtSGluZExpbWIzTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC4xMjQsXG4gICAgICAgIFwieVwiOiAtMjAuMjY3NyxcbiAgICAgICAgXCJ6XCI6IC0xLjI0MDJcbiAgICB9LFxuICAgIFwibUhpbmRMaW1iNExlZnRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0yLjQwMTYsXG4gICAgICAgIFwielwiOiA0LjQwOTRcbiAgICB9LFxuICAgIFwiZW5kX21IaW5kTGltYjRMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuMzE1LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDQuMTMzOVxuICAgIH0sXG4gICAgXCJtR3JvaW5cIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0zLjgxODksXG4gICAgICAgIFwielwiOiAyLjUxOTdcbiAgICB9LFxuICAgIFwiZW5kX21Hcm9pblwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTIuNTk4NCxcbiAgICAgICAgXCJ6XCI6IDAuMTU3NVxuICAgIH0sXG4gICAgXCJtVGFpbDFcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDEuODUwNCxcbiAgICAgICAgXCJ6XCI6IC00LjU2NjlcbiAgICB9LFxuICAgIFwibVRhaWwyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTcuNzU1OVxuICAgIH0sXG4gICAgXCJtVGFpbDNcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtNi42MTQyXG4gICAgfSxcbiAgICBcIm1UYWlsNFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC01LjU5MDVcbiAgICB9LFxuICAgIFwibVRhaWw1XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTQuNDA5NFxuICAgIH0sXG4gICAgXCJtVGFpbDZcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtMy43MDA4XG4gICAgfSxcbiAgICBcImVuZF9tVGFpbDZcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtMy41MDM5XG4gICAgfSxcbiAgICBcImxUaGlnaFwiOiB7XG4gICAgICAgIFwieFwiOiA0Ljk5MDcsXG4gICAgICAgIFwieVwiOiAtMS42MTQxLFxuICAgICAgICBcInpcIjogMS4zMjlcbiAgICB9LFxuICAgIFwibFNoaW5cIjoge1xuICAgICAgICBcInhcIjogLTEuODgzNyxcbiAgICAgICAgXCJ5XCI6IC0xOS4zMzI4LFxuICAgICAgICBcInpcIjogLTAuMDM2N1xuICAgIH0sXG4gICAgXCJsRm9vdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjA1NyxcbiAgICAgICAgXCJ5XCI6IC0yMC4yODcxLFxuICAgICAgICBcInpcIjogLTEuMTk0MVxuICAgIH0sXG4gICAgXCJtRm9vdExlZnRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0yLjM4NjYsXG4gICAgICAgIFwielwiOiA0LjQwNzdcbiAgICB9LFxuICAgIFwibVRvZUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiA0LjI5MTNcbiAgICB9LFxuICAgIFwiZW5kX21Ub2VMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC43ODc0XG4gICAgfSxcbiAgICBcInJUaGlnaFwiOiB7XG4gICAgICAgIFwieFwiOiAtNS4wNzExLFxuICAgICAgICBcInlcIjogLTEuNjE3NixcbiAgICAgICAgXCJ6XCI6IDEuMzIzNlxuICAgIH0sXG4gICAgXCJyU2hpblwiOiB7XG4gICAgICAgIFwieFwiOiAyLjAxMDUsXG4gICAgICAgIFwieVwiOiAtMTkuMzI3NixcbiAgICAgICAgXCJ6XCI6IC0wLjAzMjJcbiAgICB9LFxuICAgIFwickZvb3RcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0yMC4yODkxLFxuICAgICAgICBcInpcIjogLTEuMTkzNFxuICAgIH0sXG4gICAgXCJtRm9vdFJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMi4zODczLFxuICAgICAgICBcInpcIjogNC40MDc3XG4gICAgfSxcbiAgICBcIm1Ub2VSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDQuMjkxM1xuICAgIH0sXG4gICAgXCJlbmRfbVRvZVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC43ODc0XG4gICAgfSxcbiAgICBcIm1TcGluZTFcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDMuMzEsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcIm1TcGluZTJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0zLjMxLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJhYmRvbWVuXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAzLjMxLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJtU3BpbmUzXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiA4LjQ2OTMsXG4gICAgICAgIFwielwiOiAtMC42MDVcbiAgICB9LFxuICAgIFwibVNwaW5lNFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTguNDY5MyxcbiAgICAgICAgXCJ6XCI6IDAuNjA1XG4gICAgfSxcbiAgICBcImNoZXN0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiA4LjQ2OTMsXG4gICAgICAgIFwielwiOiAtMC42MDVcbiAgICB9LFxuICAgIFwibVdpbmdzUm9vdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC0wLjU3MzJcbiAgICB9LFxuICAgIFwibVdpbmcxUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTMuNzQwMixcbiAgICAgICAgXCJ5XCI6IDcuMTI2LFxuICAgICAgICBcInpcIjogLTUuODY2MVxuICAgIH0sXG4gICAgXCJtV2luZzJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtNi42NTM1LFxuICAgICAgICBcInlcIjogMi42Mzc4LFxuICAgICAgICBcInpcIjogLTYuNjE0MlxuICAgIH0sXG4gICAgXCJtV2luZzNSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtNy4yMDQ3LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC03LjEyNlxuICAgIH0sXG4gICAgXCJtV2luZzRGYW5SaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtNi44MTEsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTYuNzMyM1xuICAgIH0sXG4gICAgXCJlbmRfbVdpbmc0RmFuUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTIuNDQwOSxcbiAgICAgICAgXCJ5XCI6IC02LjI1OTgsXG4gICAgICAgIFwielwiOiAtMi42NzcyXG4gICAgfSxcbiAgICBcIm1XaW5nNFJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC02LjgxMSxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtNi43MzIzXG4gICAgfSxcbiAgICBcImVuZF9tV2luZzRSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtNS4xOTY4LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC01Ljc0OFxuICAgIH0sXG4gICAgXCJtV2luZzFMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDMuNzQwMixcbiAgICAgICAgXCJ5XCI6IDcuMTI2LFxuICAgICAgICBcInpcIjogLTUuODY2MVxuICAgIH0sXG4gICAgXCJtV2luZzJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDYuNjUzNSxcbiAgICAgICAgXCJ5XCI6IDIuNjM3OCxcbiAgICAgICAgXCJ6XCI6IC02LjYxNDJcbiAgICB9LFxuICAgIFwibVdpbmczTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiA3LjIwNDcsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTcuMTI2XG4gICAgfSxcbiAgICBcIm1XaW5nNEZhbkxlZnRcIjoge1xuICAgICAgICBcInhcIjogNi44MTEsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTYuNzMyM1xuICAgIH0sXG4gICAgXCJlbmRfbVdpbmc0RmFuTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAyLjQ0MDksXG4gICAgICAgIFwieVwiOiAtNi4yNTk4LFxuICAgICAgICBcInpcIjogLTIuNjc3MlxuICAgIH0sXG4gICAgXCJtV2luZzRMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDYuODExLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC02LjczMjNcbiAgICB9LFxuICAgIFwiZW5kX21XaW5nNExlZnRcIjoge1xuICAgICAgICBcInhcIjogNS4xOTY4LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC01Ljc0OFxuICAgIH0sXG4gICAgXCJyQ29sbGFyXCI6IHtcbiAgICAgICAgXCJ4XCI6IC0yLjk4MjMsXG4gICAgICAgIFwieVwiOiA2LjgzNzIsXG4gICAgICAgIFwielwiOiAtMC44NTY5XG4gICAgfSxcbiAgICBcInJTaGxkclwiOiB7XG4gICAgICAgIFwieFwiOiAtNC4zNzc0LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwickZvcmVBcm1cIjoge1xuICAgICAgICBcInhcIjogLTE0LjM1MjcsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJySGFuZFwiOiB7XG4gICAgICAgIFwieFwiOiAtMTAuMzMwNyxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcIm1IYW5kVGh1bWIxUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuMDIzNixcbiAgICAgICAgXCJ5XCI6IDAuMTU3NSxcbiAgICAgICAgXCJ6XCI6IDEuMjIwNVxuICAgIH0sXG4gICAgXCJtSGFuZFRodW1iMlJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjI1OTgsXG4gICAgICAgIFwieVwiOiAtMC4wMzk0LFxuICAgICAgICBcInpcIjogMS4xMDI0XG4gICAgfSxcbiAgICBcIm1IYW5kVGh1bWIzUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuMjIwNSxcbiAgICAgICAgXCJ5XCI6IC0wLjAzOTQsXG4gICAgICAgIFwielwiOiAwLjkwNTVcbiAgICB9LFxuICAgIFwiZW5kX21IYW5kVGh1bWIzUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuOTg0MyxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjU5MDZcbiAgICB9LFxuICAgIFwibUhhbmRQaW5reTFSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMy43NDAyLFxuICAgICAgICBcInlcIjogMC4xMTgxLFxuICAgICAgICBcInpcIjogLTEuMjIwNVxuICAgIH0sXG4gICAgXCJtSGFuZFBpbmt5MlJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjk4NDIsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogLTAuOTQ0OVxuICAgIH0sXG4gICAgXCJtSGFuZFBpbmt5M1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjcwODcsXG4gICAgICAgIFwieVwiOiAtMC4xNTc1LFxuICAgICAgICBcInpcIjogLTAuNTkwNlxuICAgIH0sXG4gICAgXCJlbmRfbUhhbmRQaW5reTNSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC42Mjk5LFxuICAgICAgICBcInlcIjogLTAuMTU3NSxcbiAgICAgICAgXCJ6XCI6IC0wLjUxMThcbiAgICB9LFxuICAgIFwibUhhbmRSaW5nMVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0zLjg5NzYsXG4gICAgICAgIFwieVwiOiAwLjM1NDMsXG4gICAgICAgIFwielwiOiAtMC4zOTM3XG4gICAgfSxcbiAgICBcIm1IYW5kUmluZzJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS40OTYxLFxuICAgICAgICBcInlcIjogLTAuMzE1LFxuICAgICAgICBcInpcIjogLTAuNTExOFxuICAgIH0sXG4gICAgXCJtSGFuZFJpbmczUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuNTc0OCxcbiAgICAgICAgXCJ5XCI6IC0wLjM1NDMsXG4gICAgICAgIFwielwiOiAtMC41MTE4XG4gICAgfSxcbiAgICBcImVuZF9tSGFuZFJpbmczUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuMTAyNCxcbiAgICAgICAgXCJ5XCI6IC0wLjIzNjIsXG4gICAgICAgIFwielwiOiAtMC4zOTM3XG4gICAgfSxcbiAgICBcIm1IYW5kSW5kZXgxUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTMuODE4OSxcbiAgICAgICAgXCJ5XCI6IDAuNTkwNSxcbiAgICAgICAgXCJ6XCI6IDEuNDk2MVxuICAgIH0sXG4gICAgXCJtSGFuZEluZGV4MlJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjQxNzMsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogMC42NjkzXG4gICAgfSxcbiAgICBcIm1IYW5kSW5kZXgzUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuMjU5OCxcbiAgICAgICAgXCJ5XCI6IC0wLjIzNjIsXG4gICAgICAgIFwielwiOiAwLjU1MTJcbiAgICB9LFxuICAgIFwiZW5kX21IYW5kSW5kZXgzUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuOTg0MyxcbiAgICAgICAgXCJ5XCI6IC0wLjE1NzUsXG4gICAgICAgIFwielwiOiAwLjQzMzFcbiAgICB9LFxuICAgIFwibUhhbmRNaWRkbGUxUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTMuOTc2NCxcbiAgICAgICAgXCJ5XCI6IDAuNTkwNSxcbiAgICAgICAgXCJ6XCI6IDAuNTExOFxuICAgIH0sXG4gICAgXCJtSGFuZE1pZGRsZTJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS41NzQ4LFxuICAgICAgICBcInlcIjogLTAuMjM2MixcbiAgICAgICAgXCJ6XCI6IC0wLjAzOTRcbiAgICB9LFxuICAgIFwibUhhbmRNaWRkbGUzUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuOTI5MSxcbiAgICAgICAgXCJ5XCI6IC0wLjMxNSxcbiAgICAgICAgXCJ6XCI6IC0wLjAzOTRcbiAgICB9LFxuICAgIFwiZW5kX21IYW5kTWlkZGxlM1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjI5OTIsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogLTAuMDc4N1xuICAgIH0sXG4gICAgXCJsQ29sbGFyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDIuOTY5LFxuICAgICAgICBcInlcIjogNi44MzcyLFxuICAgICAgICBcInpcIjogLTAuODU2OVxuICAgIH0sXG4gICAgXCJsU2hsZHJcIjoge1xuICAgICAgICBcInhcIjogNC4zNTQzLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwibEZvcmVBcm1cIjoge1xuICAgICAgICBcInhcIjogMTQuMzUyNyxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcImxIYW5kXCI6IHtcbiAgICAgICAgXCJ4XCI6IDEwLjMyMjksXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJtSGFuZFRodW1iMUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMS4wMjM2LFxuICAgICAgICBcInlcIjogMC4xNTc1LFxuICAgICAgICBcInpcIjogMS4yMjA1XG4gICAgfSxcbiAgICBcIm1IYW5kVGh1bWIyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjI1OTgsXG4gICAgICAgIFwieVwiOiAtMC4wMzk0LFxuICAgICAgICBcInpcIjogMS4xMDI0XG4gICAgfSxcbiAgICBcIm1IYW5kVGh1bWIzTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjIyMDUsXG4gICAgICAgIFwieVwiOiAtMC4wMzk0LFxuICAgICAgICBcInpcIjogMC45MDU1XG4gICAgfSxcbiAgICBcImVuZF9tSGFuZFRodW1iM0xlZnRcIjoge1xuICAgICAgICBcInhcIjogMC45ODQzLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuNTkwNlxuICAgIH0sXG4gICAgXCJtSGFuZFBpbmt5MUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMy43NDAyLFxuICAgICAgICBcInlcIjogMC4xMTgxLFxuICAgICAgICBcInpcIjogLTEuMjIwNVxuICAgIH0sXG4gICAgXCJtSGFuZFBpbmt5MkxlZnRcIjoge1xuICAgICAgICBcInhcIjogMC45ODQyLFxuICAgICAgICBcInlcIjogLTAuMjM2MixcbiAgICAgICAgXCJ6XCI6IC0wLjk0NDlcbiAgICB9LFxuICAgIFwibUhhbmRQaW5reTNMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuNzA4NyxcbiAgICAgICAgXCJ5XCI6IC0wLjE1NzUsXG4gICAgICAgIFwielwiOiAtMC41OTA2XG4gICAgfSxcbiAgICBcImVuZF9tSGFuZFBpbmt5M0xlZnRcIjoge1xuICAgICAgICBcInhcIjogMC42Mjk5LFxuICAgICAgICBcInlcIjogLTAuMTU3NSxcbiAgICAgICAgXCJ6XCI6IC0wLjUxMThcbiAgICB9LFxuICAgIFwibUhhbmRSaW5nMUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMy44OTc2LFxuICAgICAgICBcInlcIjogMC4zNTQzLFxuICAgICAgICBcInpcIjogLTAuMzkzN1xuICAgIH0sXG4gICAgXCJtSGFuZFJpbmcyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjQ5NjEsXG4gICAgICAgIFwieVwiOiAtMC4zMTUsXG4gICAgICAgIFwielwiOiAtMC41MTE4XG4gICAgfSxcbiAgICBcIm1IYW5kUmluZzNMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuNTc0OCxcbiAgICAgICAgXCJ5XCI6IC0wLjM1NDMsXG4gICAgICAgIFwielwiOiAtMC41MTE4XG4gICAgfSxcbiAgICBcImVuZF9tSGFuZFJpbmczTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjEwMjQsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogLTAuMzkzN1xuICAgIH0sXG4gICAgXCJtSGFuZEluZGV4MUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMy44MTg5LFxuICAgICAgICBcInlcIjogMC41OTA1LFxuICAgICAgICBcInpcIjogMS40OTYxXG4gICAgfSxcbiAgICBcIm1IYW5kSW5kZXgyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjQxNzMsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogMC42NjkzXG4gICAgfSxcbiAgICBcIm1IYW5kSW5kZXgzTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjI1OTgsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogMC41NTEyXG4gICAgfSxcbiAgICBcImVuZF9tSGFuZEluZGV4M0xlZnRcIjoge1xuICAgICAgICBcInhcIjogMC45ODQzLFxuICAgICAgICBcInlcIjogLTAuMTU3NSxcbiAgICAgICAgXCJ6XCI6IDAuNDMzMVxuICAgIH0sXG4gICAgXCJtSGFuZE1pZGRsZTFMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDMuOTc2NCxcbiAgICAgICAgXCJ5XCI6IDAuNTkwNSxcbiAgICAgICAgXCJ6XCI6IDAuNTExOFxuICAgIH0sXG4gICAgXCJtSGFuZE1pZGRsZTJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuNTc0OCxcbiAgICAgICAgXCJ5XCI6IC0wLjIzNjIsXG4gICAgICAgIFwielwiOiAtMC4wMzk0XG4gICAgfSxcbiAgICBcIm1IYW5kTWlkZGxlM0xlZnRcIjoge1xuICAgICAgICBcInhcIjogMS45MjkxLFxuICAgICAgICBcInlcIjogLTAuMzE1LFxuICAgICAgICBcInpcIjogLTAuMDM5NFxuICAgIH0sXG4gICAgXCJlbmRfbUhhbmRNaWRkbGUzTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjI5OTIsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogLTAuMDc4N1xuICAgIH0sXG4gICAgXCJuZWNrXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAxMC4zODA0LFxuICAgICAgICBcInpcIjogLTAuMzg5M1xuICAgIH0sXG4gICAgXCJoZWFkXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAzLjU3MzEsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcIm1GYWNlUm9vdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMS42OTM0LFxuICAgICAgICBcInpcIjogMC45MTA0XG4gICAgfSxcbiAgICBcIm1GYWNlTm9zZUJyaWRnZVwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMC43MTE4LFxuICAgICAgICBcInpcIjogMy4zMTRcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlTm9zZUJyaWRnZVwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMC4zMTUsXG4gICAgICAgIFwielwiOiAwLjU5MDZcbiAgICB9LFxuICAgIFwibUZhY2VFeWVjb3JuZXJJbm5lclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjcxMzUsXG4gICAgICAgIFwieVwiOiAxLjExMjIsXG4gICAgICAgIFwielwiOiAyLjczMTNcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRXllY29ybmVySW5uZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuNjI5OVxuICAgIH0sXG4gICAgXCJtRmFjZUV5ZWNvcm5lcklubmVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjcxMzUsXG4gICAgICAgIFwieVwiOiAxLjExMjIsXG4gICAgICAgIFwielwiOiAyLjczMTNcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRXllY29ybmVySW5uZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC42Mjk5XG4gICAgfSxcbiAgICBcIm1GYWNlVGVldGhVcHBlclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTEuMTQ3OCxcbiAgICAgICAgXCJ6XCI6IDAuNzI4M1xuICAgIH0sXG4gICAgXCJtRmFjZUxpcFVwcGVyQ2VudGVyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMC4wOTY0LFxuICAgICAgICBcInpcIjogMS42Mzg4XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUxpcFVwcGVyQ2VudGVyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLjA3ODcsXG4gICAgICAgIFwielwiOiAxLjY5MjlcbiAgICB9LFxuICAgIFwibUZhY2VMaXBDb3JuZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjY5MTksXG4gICAgICAgIFwieVwiOiAtMC4zNjQyLFxuICAgICAgICBcInpcIjogMS4wMTk3XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUxpcENvcm5lclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0yLjAwNzksXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMS43NzE3XG4gICAgfSxcbiAgICBcIm1GYWNlTGlwQ29ybmVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC42OTE5LFxuICAgICAgICBcInlcIjogLTAuMzY0MixcbiAgICAgICAgXCJ6XCI6IDEuMDE5N1xuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VMaXBDb3JuZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDIuMDA3OSxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAxLjc3MTdcbiAgICB9LFxuICAgIFwibUZhY2VMaXBVcHBlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMC4xNDc4LFxuICAgICAgICBcInpcIjogMS42Mzg4XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUxpcFVwcGVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuNTkwNixcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAxLjYxNDJcbiAgICB9LFxuICAgIFwibUZhY2VMaXBVcHBlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0wLjE0NzgsXG4gICAgICAgIFwielwiOiAxLjYzODhcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlTGlwVXBwZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuNTkwNixcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAxLjYxNDJcbiAgICB9LFxuICAgIFwibUZhY2VOb3NlQmFzZVwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTAuNzA0NyxcbiAgICAgICAgXCJ6XCI6IDMuNDIzMlxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VOb3NlQmFzZVwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuNTUxMlxuICAgIH0sXG4gICAgXCJtRmFjZUZvcmVoZWFkQ2VudGVyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAyLjMwMjcsXG4gICAgICAgIFwielwiOiAyLjUyMzdcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRm9yZWhlYWRDZW50ZXJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAxLjQxNzNcbiAgICB9LFxuICAgIFwibUZhY2VKYXdTaGFwZXJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUphd1NoYXBlclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC0wLjY2OTNcbiAgICB9LFxuICAgIFwibUZhY2VKYXdcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0wLjUzMzksXG4gICAgICAgIFwielwiOiAtMC4wMzY0XG4gICAgfSxcbiAgICBcIm1GYWNlVGVldGhMb3dlclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTEuNDg2OSxcbiAgICAgICAgXCJ6XCI6IDAuNzY0OFxuICAgIH0sXG4gICAgXCJtRmFjZVRvbmd1ZUJhc2VcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAuMTgyMSxcbiAgICAgICAgXCJ6XCI6IDEuNDIwM1xuICAgIH0sXG4gICAgXCJtRmFjZVRvbmd1ZVRpcFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMC4yNTQ5LFxuICAgICAgICBcInpcIjogMC44MDEyXG4gICAgfSxcbiAgICBcImVuZF9tRmFjZVRvbmd1ZVRpcFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuMzkzN1xuICAgIH0sXG4gICAgXCJtRmFjZUxpcExvd2VyQ2VudGVyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMS42Mzg4XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUxpcExvd2VyQ2VudGVyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLjA3ODcsXG4gICAgICAgIFwielwiOiAxLjU3NDhcbiAgICB9LFxuICAgIFwibUZhY2VMaXBMb3dlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMS42Mzg4XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUxpcExvd2VyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuNjY5MyxcbiAgICAgICAgXCJ5XCI6IDAuMTk2OSxcbiAgICAgICAgXCJ6XCI6IDEuMzM4NlxuICAgIH0sXG4gICAgXCJtRmFjZUxpcExvd2VyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDEuNjM4OFxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VMaXBMb3dlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMC42NjkzLFxuICAgICAgICBcInlcIjogMC4xOTY5LFxuICAgICAgICBcInpcIjogMS4zMzg2XG4gICAgfSxcbiAgICBcIm1GYWNlQ2hpblwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTIuMDIxMixcbiAgICAgICAgXCJ6XCI6IDIuNTMxXG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUNoaW5cIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0wLjcwODcsXG4gICAgICAgIFwielwiOiAwLjgyNjhcbiAgICB9LFxuICAgIFwibUZhY2VDaGVla1VwcGVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuMjY2MyxcbiAgICAgICAgXCJ5XCI6IC0wLjIzMTMsXG4gICAgICAgIFwielwiOiAyLjU0OTJcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlQ2hlZWtVcHBlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjU5MDYsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC44NjYxXG4gICAgfSxcbiAgICBcIm1GYWNlQ2hlZWtMb3dlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjI2NjMsXG4gICAgICAgIFwieVwiOiAtMS4xOTQxLFxuICAgICAgICBcInpcIjogMS44MjA5XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUNoZWVrTG93ZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS4xODExLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuNTExOFxuICAgIH0sXG4gICAgXCJtRmFjZUNoZWVrVXBwZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuMjY2MyxcbiAgICAgICAgXCJ5XCI6IC0wLjIzMTMsXG4gICAgICAgIFwielwiOiAyLjU0OTJcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlQ2hlZWtVcHBlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMC41OTA2LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuODY2MVxuICAgIH0sXG4gICAgXCJtRmFjZUNoZWVrTG93ZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuMjY2MyxcbiAgICAgICAgXCJ5XCI6IC0xLjE5NDEsXG4gICAgICAgIFwielwiOiAxLjgyMDlcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlQ2hlZWtMb3dlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMS4xODExLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuNTExOFxuICAgIH0sXG4gICAgXCJtRmFjZU5vc2VSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC41NTg3LFxuICAgICAgICBcInlcIjogLTAuMTk1NyxcbiAgICAgICAgXCJ6XCI6IDMuMTMxOVxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VOb3NlUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuMTU3NSxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjU5MDZcbiAgICB9LFxuICAgIFwibUZhY2VOb3NlQ2VudGVyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMC4wNTM0LFxuICAgICAgICBcInpcIjogMy43MTQ2XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZU5vc2VDZW50ZXJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjk4NDNcbiAgICB9LFxuICAgIFwibUZhY2VOb3NlTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjU1ODcsXG4gICAgICAgIFwieVwiOiAtMC4xOTU3LFxuICAgICAgICBcInpcIjogMy4xMzE5XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZU5vc2VMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuMTU3NSxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjU5MDZcbiAgICB9LFxuICAgIFwibUZhY2VFYXIxUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTIuOTc5NSxcbiAgICAgICAgXCJ5XCI6IDAuMDcxMixcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwibUZhY2VFYXIyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuNzA4NyxcbiAgICAgICAgXCJ5XCI6IDAuOTg0MixcbiAgICAgICAgXCJ6XCI6IC0wLjc0OFxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFYXIyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDEuMjk5MixcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwibUZhY2VFYXIxTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAyLjk3OTUsXG4gICAgICAgIFwieVwiOiAwLjA3MTIsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcIm1GYWNlRWFyMkxlZnRcIjoge1xuICAgICAgICBcInhcIjogMC43MDg3LFxuICAgICAgICBcInlcIjogMC45ODQyLFxuICAgICAgICBcInpcIjogLTAuNzQ4XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUVhcjJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAxLjI5OTIsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcIm1GYWNlRXllTGlkTG93ZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS40MTU2LFxuICAgICAgICBcInlcIjogMS4xODM0LFxuICAgICAgICBcInpcIjogMi42NTg1XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUV5ZUxpZExvd2VyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0wLjI3NTYsXG4gICAgICAgIFwielwiOiAwLjk0NDlcbiAgICB9LFxuICAgIFwibUZhY2VFeWVMaWRVcHBlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjQxNTYsXG4gICAgICAgIFwieVwiOiAxLjE4ODcsXG4gICAgICAgIFwielwiOiAyLjY1ODVcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRXllTGlkVXBwZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMC4xOTY5LFxuICAgICAgICBcInpcIjogMS4wNjNcbiAgICB9LFxuICAgIFwibUZhY2VFeWVMaWRMb3dlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMS40MTU2LFxuICAgICAgICBcInlcIjogMS4xODM0LFxuICAgICAgICBcInpcIjogMi42NTg1XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUV5ZUxpZExvd2VyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTAuMjc1NixcbiAgICAgICAgXCJ6XCI6IDAuOTQ0OVxuICAgIH0sXG4gICAgXCJtRmFjZUV5ZUxpZFVwcGVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjQxNTYsXG4gICAgICAgIFwieVwiOiAxLjE4ODcsXG4gICAgICAgIFwielwiOiAyLjY1ODVcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRXllTGlkVXBwZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLjE5NjksXG4gICAgICAgIFwielwiOiAxLjA2M1xuICAgIH0sXG4gICAgXCJtRmFjZUV5ZWJyb3dJbm5lclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjkwNixcbiAgICAgICAgXCJ5XCI6IDEuODU3OCxcbiAgICAgICAgXCJ6XCI6IDIuNzMxM1xuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVicm93SW5uZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDEuMDIzNlxuICAgIH0sXG4gICAgXCJtRmFjZUV5ZWJyb3dDZW50ZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS42NzExLFxuICAgICAgICBcInlcIjogMi4xNDY3LFxuICAgICAgICBcInpcIjogMi41NDkyXG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUV5ZWJyb3dDZW50ZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDEuMDYzXG4gICAgfSxcbiAgICBcIm1GYWNlRXllYnJvd091dGVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuOTE2OCxcbiAgICAgICAgXCJ5XCI6IDEuNzI5LFxuICAgICAgICBcInpcIjogMi4zMzA3XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUV5ZWJyb3dPdXRlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjUxMTgsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC45MDU1XG4gICAgfSxcbiAgICBcIm1GYWNlRXllYnJvd0lubmVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjkwNixcbiAgICAgICAgXCJ5XCI6IDEuODU3OCxcbiAgICAgICAgXCJ6XCI6IDIuNzMxM1xuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVicm93SW5uZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMS4wMjM2XG4gICAgfSxcbiAgICBcIm1GYWNlRXllYnJvd0NlbnRlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMS42NzExLFxuICAgICAgICBcInlcIjogMi4xNDY3LFxuICAgICAgICBcInpcIjogMi41NDkyXG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUV5ZWJyb3dDZW50ZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMS4wNjNcbiAgICB9LFxuICAgIFwibUZhY2VFeWVicm93T3V0ZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuOTE2OCxcbiAgICAgICAgXCJ5XCI6IDEuNzI5LFxuICAgICAgICBcInpcIjogMi4zMzA3XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUV5ZWJyb3dPdXRlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMC41MTE4LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuOTA1NVxuICAgIH0sXG4gICAgXCJtRmFjZUZvcmVoZWFkUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuMzAzNSxcbiAgICAgICAgXCJ5XCI6IDMuMDYwOCxcbiAgICAgICAgXCJ6XCI6IDIuMzMwN1xuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VGb3JlaGVhZFJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjE1NzUsXG4gICAgICAgIFwieVwiOiAwLjcwODcsXG4gICAgICAgIFwielwiOiAwLjk0NDlcbiAgICB9LFxuICAgIFwibUZhY2VGb3JlaGVhZExlZnRcIjoge1xuICAgICAgICBcInhcIjogMS4zMDM1LFxuICAgICAgICBcInlcIjogMy4wNjA4LFxuICAgICAgICBcInpcIjogMi4zMzA3XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUZvcmVoZWFkTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjE1NzUsXG4gICAgICAgIFwieVwiOiAwLjcwODcsXG4gICAgICAgIFwielwiOiAwLjk0NDlcbiAgICB9LFxuICAgIFwibUZhY2VFeWVBbHRMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuNDE1NixcbiAgICAgICAgXCJ5XCI6IDEuMTgzNixcbiAgICAgICAgXCJ6XCI6IDIuNjc1MlxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVBbHRMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC45ODQzXG4gICAgfSxcbiAgICBcIm1GYWNlRXllQWx0UmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuNDE1NixcbiAgICAgICAgXCJ5XCI6IDEuMTgzNixcbiAgICAgICAgXCJ6XCI6IDIuNjc1NFxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVBbHRSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuOTg0M1xuICAgIH0sXG4gICAgXCJtRXllTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjQyMDMsXG4gICAgICAgIFwieVwiOiAyLjg3NyxcbiAgICAgICAgXCJ6XCI6IDMuNTg1N1xuICAgIH0sXG4gICAgXCJlbmRfbUV5ZUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjk4NDNcbiAgICB9LFxuICAgIFwibUV5ZVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjQyMDMsXG4gICAgICAgIFwieVwiOiAyLjg3NyxcbiAgICAgICAgXCJ6XCI6IDMuNTg1OVxuICAgIH0sXG4gICAgXCJlbmRfbUV5ZVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC45ODQzXG4gICAgfSxcbiAgICBcImZpZ3VyZUhhaXJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDIuNjAzOCxcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwiZW5kX2ZpZ3VyZUhhaXJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDEuMjk5MixcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9XG59O1xuIiwiXCJ1c2Ugc3RyaWN0XCI7XG5PYmplY3QuZGVmaW5lUHJvcGVydHkoZXhwb3J0cywgXCJfX2VzTW9kdWxlXCIsIHsgdmFsdWU6IHRydWUgfSk7XG5leHBvcnRzLnBhcnNlQW5pbSA9IHBhcnNlQW5pbTtcbmV4cG9ydHMuZGlzdHJpYnV0ZVNpbmdsZUZyYW1lID0gZGlzdHJpYnV0ZVNpbmdsZUZyYW1lO1xuZXhwb3J0cy5wYXJzZUJWSCA9IHBhcnNlQlZIO1xuY29uc3QgdXRpbHNfMSA9IHJlcXVpcmUoXCIuL3V0aWxzXCIpO1xuZnVuY3Rpb24gcGFyc2VBbmltKGFycmF5QnVmZmVyKSB7XG4gICAgY29uc3QgdmlldyA9IG5ldyBEYXRhVmlldyhhcnJheUJ1ZmZlcik7XG4gICAgbGV0IG9mZnNldCA9IDA7XG4gICAgZnVuY3Rpb24gcmVhZFUxNigpIHtcbiAgICAgICAgY29uc3QgdmFsdWUgPSB2aWV3LmdldFVpbnQxNihvZmZzZXQsIHRydWUpO1xuICAgICAgICBvZmZzZXQgKz0gMjtcbiAgICAgICAgcmV0dXJuIHZhbHVlO1xuICAgIH1cbiAgICBmdW5jdGlvbiByZWFkVTMyKCkge1xuICAgICAgICBjb25zdCB2YWx1ZSA9IHZpZXcuZ2V0VWludDE2KG9mZnNldCwgdHJ1ZSk7XG4gICAgICAgIG9mZnNldCArPSA0O1xuICAgICAgICByZXR1cm4gdmFsdWU7XG4gICAgfVxuICAgIGZ1bmN0aW9uIHJlYWRTMzIoKSB7XG4gICAgICAgIGNvbnN0IHZhbHVlID0gdmlldy5nZXRJbnQzMihvZmZzZXQsIHRydWUpO1xuICAgICAgICBvZmZzZXQgKz0gNDtcbiAgICAgICAgcmV0dXJuIHZhbHVlO1xuICAgIH1cbiAgICBmdW5jdGlvbiByZWFkRjMyKCkge1xuICAgICAgICBjb25zdCB2YWx1ZSA9IHZpZXcuZ2V0RmxvYXQzMihvZmZzZXQsIHRydWUpO1xuICAgICAgICBvZmZzZXQgKz0gNDtcbiAgICAgICAgcmV0dXJuIHZhbHVlO1xuICAgIH1cbiAgICBmdW5jdGlvbiByZWFkU3RyaW5nKCkge1xuICAgICAgICBsZXQgc3RyID0gXCJcIjtcbiAgICAgICAgd2hpbGUgKG9mZnNldCA8IHZpZXcuYnl0ZUxlbmd0aCkge1xuICAgICAgICAgICAgY29uc3QgY2hhciA9IHZpZXcuZ2V0VWludDgob2Zmc2V0KyspO1xuICAgICAgICAgICAgaWYgKGNoYXIgPT09IDApXG4gICAgICAgICAgICAgICAgYnJlYWs7IC8vIE5VTEwtdGVybWluYXRlZFxuICAgICAgICAgICAgc3RyICs9IFN0cmluZy5mcm9tQ2hhckNvZGUoY2hhcik7XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIHN0cjtcbiAgICB9XG4gICAgZnVuY3Rpb24gcmVhZEZpeGVkU3RyaW5nKHNpemUpIHtcbiAgICAgICAgbGV0IHN0ciA9IFwiXCI7XG4gICAgICAgIGZvciAobGV0IGkgPSAwOyBpIDwgc2l6ZTsgaSsrKSB7XG4gICAgICAgICAgICBjb25zdCBjaGFyID0gdmlldy5nZXRVaW50OChvZmZzZXQrKyk7XG4gICAgICAgICAgICBpZiAoY2hhciAhPT0gMClcbiAgICAgICAgICAgICAgICBzdHIgKz0gU3RyaW5nLmZyb21DaGFyQ29kZShjaGFyKTtcbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gc3RyO1xuICAgIH1cbiAgICBjb25zdCB2ZXJzaW9uID0gcmVhZFUxNigpO1xuICAgIGNvbnN0IHN1Yl92ZXJzaW9uID0gcmVhZFUxNigpO1xuICAgIGNvbnN0IGJhc2VfcHJpb3JpdHkgPSByZWFkUzMyKCk7XG4gICAgY29uc3QgZHVyYXRpb24gPSByZWFkRjMyKCk7XG4gICAgY29uc3QgZW1vdGVfbmFtZSA9IHJlYWRTdHJpbmcoKTtcbiAgICBjb25zdCBsb29wX2luX3BvaW50ID0gcmVhZEYzMigpO1xuICAgIGNvbnN0IGxvb3Bfb3V0X3BvaW50ID0gcmVhZEYzMigpO1xuICAgIGNvbnN0IGxvb3AgPSByZWFkUzMyKCk7XG4gICAgY29uc3QgZWFzZV9pbl9kdXJhdGlvbiA9IHJlYWRGMzIoKTtcbiAgICBjb25zdCBlYXNlX291dF9kdXJhdGlvbiA9IHJlYWRGMzIoKTtcbiAgICBjb25zdCBoYW5kX3Bvc2UgPSByZWFkVTMyKCk7XG4gICAgY29uc3QgbnVtX2pvaW50cyA9IHJlYWRVMzIoKTtcbiAgICBjb25zdCBqb2ludHMgPSBbXTtcbiAgICBmb3IgKGxldCBpID0gMDsgaSA8IG51bV9qb2ludHM7IGkrKykge1xuICAgICAgICBjb25zdCBqb2ludF9uYW1lID0gcmVhZFN0cmluZygpO1xuICAgICAgICBjb25zdCBqb2ludF9wcmlvcml0eSA9IHJlYWRTMzIoKTtcbiAgICAgICAgY29uc3QgbnVtX3JvdF9rZXlzID0gcmVhZFMzMigpO1xuICAgICAgICBjb25zdCByb3RhdGlvbl9rZXlzID0gW107XG4gICAgICAgIGZvciAobGV0IGogPSAwOyBqIDwgbnVtX3JvdF9rZXlzOyBqKyspIHtcbiAgICAgICAgICAgIGNvbnN0IHRpbWUgPSByZWFkVTE2KCk7XG4gICAgICAgICAgICBjb25zdCByb3RfeCA9IHJlYWRVMTYoKTtcbiAgICAgICAgICAgIGNvbnN0IHJvdF95ID0gcmVhZFUxNigpO1xuICAgICAgICAgICAgY29uc3Qgcm90X3ogPSByZWFkVTE2KCk7XG4gICAgICAgICAgICByb3RhdGlvbl9rZXlzLnB1c2goeyB0aW1lOiBpbnRUb0Zsb2F0KHRpbWUsIDAsIGR1cmF0aW9uKSwgeDogaW50VG9GbG9hdChyb3RfeCwgLTEsIDEpLCB5OiBpbnRUb0Zsb2F0KHJvdF95LCAtMSwgMSksIHo6IGludFRvRmxvYXQocm90X3osIC0xLCAxKSB9KTtcbiAgICAgICAgfVxuICAgICAgICBjb25zdCBudW1fcG9zX2tleXMgPSByZWFkUzMyKCk7XG4gICAgICAgIGNvbnN0IHBvc2l0aW9uX2tleXMgPSBbXTtcbiAgICAgICAgZm9yIChsZXQgaiA9IDA7IGogPCBudW1fcG9zX2tleXM7IGorKykge1xuICAgICAgICAgICAgY29uc3QgdGltZSA9IHJlYWRVMTYoKTtcbiAgICAgICAgICAgIGNvbnN0IHBvc194ID0gcmVhZFUxNigpO1xuICAgICAgICAgICAgY29uc3QgcG9zX3kgPSByZWFkVTE2KCk7XG4gICAgICAgICAgICBjb25zdCBwb3NfeiA9IHJlYWRVMTYoKTtcbiAgICAgICAgICAgIHBvc2l0aW9uX2tleXMucHVzaCh7IHRpbWU6IGludFRvRmxvYXQodGltZSwgMCwgZHVyYXRpb24pLCB4OiBpbnRUb0Zsb2F0KHBvc194LCAtNSwgNSksIHk6IGludFRvRmxvYXQocG9zX3ksIC01LCA1KSwgejogaW50VG9GbG9hdChwb3NfeiwgLTUsIDUpIH0pO1xuICAgICAgICB9XG4gICAgICAgIGpvaW50cy5wdXNoKHsgam9pbnRfbmFtZSwgam9pbnRfcHJpb3JpdHksIHJvdGF0aW9uX2tleXMsIHBvc2l0aW9uX2tleXMgfSk7XG4gICAgfVxuICAgIGNvbnN0IG51bV9jb25zdHJhaW50cyA9IHJlYWRTMzIoKTtcbiAgICBjb25zdCBjb25zdHJhaW50cyA9IFtdO1xuICAgIGZvciAobGV0IGkgPSAwOyBpIDwgbnVtX2NvbnN0cmFpbnRzOyBpKyspIHtcbiAgICAgICAgY29uc3QgY2hhaW5fbGVuZ3RoID0gdmlldy5nZXRVaW50OChvZmZzZXQrKyk7XG4gICAgICAgIGNvbnN0IGNvbnN0cmFpbnRfdHlwZSA9IHZpZXcuZ2V0VWludDgob2Zmc2V0KyspO1xuICAgICAgICBjb25zdCBzb3VyY2Vfdm9sdW1lID0gcmVhZEZpeGVkU3RyaW5nKDE2KTtcbiAgICAgICAgY29uc3Qgc291cmNlX29mZnNldCA9IFtyZWFkRjMyKCksIHJlYWRGMzIoKSwgcmVhZEYzMigpXTtcbiAgICAgICAgY29uc3QgdGFyZ2V0X3ZvbHVtZSA9IHJlYWRGaXhlZFN0cmluZygxNik7XG4gICAgICAgIGNvbnN0IHRhcmdldF9vZmZzZXQgPSBbcmVhZEYzMigpLCByZWFkRjMyKCksIHJlYWRGMzIoKV07XG4gICAgICAgIGNvbnN0IHRhcmdldF9kaXIgPSBbcmVhZEYzMigpLCByZWFkRjMyKCksIHJlYWRGMzIoKV07XG4gICAgICAgIGNvbnN0IGVhc2VfaW5fc3RhcnQgPSByZWFkRjMyKCk7XG4gICAgICAgIGNvbnN0IGVhc2VfaW5fc3RvcCA9IHJlYWRGMzIoKTtcbiAgICAgICAgY29uc3QgZWFzZV9vdXRfc3RhcnQgPSByZWFkRjMyKCk7XG4gICAgICAgIGNvbnN0IGVhc2Vfb3V0X3N0b3AgPSByZWFkRjMyKCk7XG4gICAgICAgIGNvbnN0cmFpbnRzLnB1c2goe1xuICAgICAgICAgICAgY2hhaW5fbGVuZ3RoLCBjb25zdHJhaW50X3R5cGUsIHNvdXJjZV92b2x1bWUsIHNvdXJjZV9vZmZzZXQsXG4gICAgICAgICAgICB0YXJnZXRfdm9sdW1lLCB0YXJnZXRfb2Zmc2V0LCB0YXJnZXRfZGlyLCBlYXNlX2luX3N0YXJ0LCBlYXNlX2luX3N0b3AsXG4gICAgICAgICAgICBlYXNlX291dF9zdGFydCwgZWFzZV9vdXRfc3RvcFxuICAgICAgICB9KTtcbiAgICB9XG4gICAgam9pbnRzLmZvckVhY2goKGl0ZW0pID0+IGl0ZW0ucm90YXRpb25fa2V5cy5mb3JFYWNoKChyb3QpID0+IHtcbiAgICAgICAgaWYgKCFpdGVtLmV1bGVyX2tleXMpIHtcbiAgICAgICAgICAgIGl0ZW0uZXVsZXJfa2V5cyA9IFtdO1xuICAgICAgICB9XG4gICAgICAgIGl0ZW0uZXVsZXJfa2V5cy5wdXNoKCgwLCB1dGlsc18xLnRvRXVsZXJzKShyb3QpKTtcbiAgICB9KSk7XG4gICAgcmV0dXJuIHsgdmVyc2lvbiwgc3ViX3ZlcnNpb24sIGR1cmF0aW9uLCBlbW90ZV9uYW1lLCBsb29wLCBqb2ludHMsIGNvbnN0cmFpbnRzIH07XG59XG5mdW5jdGlvbiBpbnRUb0Zsb2F0KHZhbCwgbWluLCBtYXgpIHtcbiAgICBjb25zdCBvbmUgPSAobWF4IC0gbWluKSAvIDY1NTM1LjA7XG4gICAgY29uc3QgcmVzdWx0ID0gbWluICsgdmFsICogb25lO1xuICAgIGlmIChNYXRoLmFicyhyZXN1bHQpIDwgb25lKSB7XG4gICAgICAgIHJldHVybiAwO1xuICAgIH1cbiAgICByZXR1cm4gcmVzdWx0O1xufVxuZnVuY3Rpb24gZW51bWVyYXRlKGNvbnRlbnQsIGtleSwgYWx0ZXIpIHtcbiAgICBsZXQgcmVzdWx0ID0gY29udGVudDtcbiAgICBsZXQgY291bnQgPSAwO1xuICAgIHdoaWxlIChyZXN1bHQuaW5jbHVkZXMoXCJcXFwiXCIgKyBrZXkgKyBcIlxcXCJcIikpIHtcbiAgICAgICAgcmVzdWx0ID0gcmVzdWx0LnJlcGxhY2UoXCJcXFwiXCIgKyBrZXkgKyBcIlxcXCJcIiwgXCJcXFwiXCIgKyBhbHRlciArIGNvdW50ICsgXCJcXFwiXCIpO1xuICAgICAgICBjb3VudCArPSAxO1xuICAgIH1cbiAgICByZXR1cm4gcmVzdWx0O1xufVxuZnVuY3Rpb24gY29tcG9zZShub2RlLCBuYW1lKSB7XG4gICAgY29uc3QgY2hpbGRyZW4gPSBbXTtcbiAgICBjb25zdCByZXN1bHQgPSB7fTtcbiAgICBsZXQgY250cyA9IFtdO1xuICAgIGxldCBqbnRzID0gW107XG4gICAgaWYgKE9iamVjdC5rZXlzKG5vZGUpLmluY2x1ZGVzKFwiRW5kIFNpdGVcIikpIHtcbiAgICAgICAgbm9kZS5qbnQxID0gXCJlbmRcIjtcbiAgICB9XG4gICAgT2JqZWN0LmtleXMobm9kZSkuZm9yRWFjaChpdGVtID0+IHtcbiAgICAgICAgaWYgKGl0ZW0uaW5jbHVkZXMoXCJjbnRcIikpIHtcbiAgICAgICAgICAgIGNudHMucHVzaChpdGVtKTtcbiAgICAgICAgICAgIHJldHVybjtcbiAgICAgICAgfVxuICAgICAgICBpZiAoaXRlbS5pbmNsdWRlcyhcImpudFwiKSkge1xuICAgICAgICAgICAgam50cy5wdXNoKGl0ZW0pO1xuICAgICAgICB9XG4gICAgfSk7XG4gICAgY250cyA9IGNudHMuc29ydCgodmFsMSwgdmFsMikgPT4geyByZXR1cm4gcGFyc2VJbnQodmFsMS5yZXBsYWNlKFwiY250XCIsIFwiXCIpKSAtIHBhcnNlSW50KHZhbDIucmVwbGFjZShcImNudFwiLCBcIlwiKSk7IH0pO1xuICAgIGpudHMgPSBqbnRzLnNvcnQoKHZhbDEsIHZhbDIpID0+IHsgcmV0dXJuIHBhcnNlSW50KHZhbDEucmVwbGFjZShcImpudFwiLCBcIlwiKSkgLSBwYXJzZUludCh2YWwyLnJlcGxhY2UoXCJqbnRcIiwgXCJcIikpOyB9KTtcbiAgICBqbnRzLmZvckVhY2goKGl0ZW0sIGkpID0+IGNoaWxkcmVuLnB1c2goY29tcG9zZShub2RlW2NudHNbaV1dLCBub2RlW2l0ZW1dLnRyaW0oKSkpKTtcbiAgICBpZiAobm9kZS5PRkZTRVQpIHtcbiAgICAgICAgY29uc3Qgb2Zmc2V0ID0gbm9kZS5PRkZTRVQudHJpbSgpLnNwbGl0KFwiIFwiKS5maWx0ZXIoKGl0ZW0pID0+ICFpc05hTihwYXJzZUZsb2F0KGl0ZW0pKSkubWFwKHBhcnNlRmxvYXQpO1xuICAgICAgICByZXN1bHQub2Zmc2V0ID0geyB4OiBvZmZzZXRbMF0sIHk6IG9mZnNldFsxXSwgejogb2Zmc2V0WzJdIH07XG4gICAgfVxuICAgIGlmIChub2RlLkNIQU5ORUxTKSB7XG4gICAgICAgIGNvbnN0IGNoYW5uZWxzID0gbm9kZS5DSEFOTkVMUy50cmltKCkuc3BsaXQoXCIgXCIpLm1hcCgoaXRlbSkgPT4gaXRlbS50cmltKCkpLmZpbHRlcigoaXRlbSkgPT4gISFpdGVtKS5maWx0ZXIoKGl0ZW0pID0+IGlzTmFOKHBhcnNlRmxvYXQoaXRlbSkpKTtcbiAgICAgICAgcmVzdWx0LmNoYW5uZWxzID0gY2hhbm5lbHM7XG4gICAgfVxuICAgIHJlc3VsdC5idmhOYW1lID0gbmFtZTtcbiAgICBpZiAoY2hpbGRyZW4ubGVuZ3RoKSB7XG4gICAgICAgIHJlc3VsdC5jaGlsZHJlbiA9IGNoaWxkcmVuO1xuICAgIH1cbiAgICByZXR1cm4gcmVzdWx0O1xufVxuZnVuY3Rpb24gY2xlYW51cChkYXRhKSB7XG4gICAgZGF0YS5qbnQwID0gZGF0YS5ST09UO1xuICAgIHJldHVybiBjb21wb3NlKGRhdGEsIFwicm9vdFwiKTtcbn1cbmZ1bmN0aW9uIHBhcnNlRnJhbWVzKHJvd3MpIHtcbiAgICBjb25zdCBzcGxpdGVkUm93cyA9IHJvd3MubWFwKGl0ZW0gPT4gaXRlbS5zcGxpdChcIiBcIikubWFwKGl0ZW0gPT4gaXRlbS50cmltKCkpLmZpbHRlcihpdGVtID0+ICEhaXRlbSkpO1xuICAgIHJldHVybiBzcGxpdGVkUm93cy5tYXAoaXRlbSA9PiBpdGVtLm1hcChwYXJzZUZsb2F0KSk7XG59XG5mdW5jdGlvbiBwYXJzZUZyYW1lc1BhcnQoZnJhbWVzUGFydCkge1xuICAgIGNvbnN0IGZyYW1lc1Jvd3MgPSBmcmFtZXNQYXJ0LnNwbGl0KFwiXFxuXCIpO1xuICAgIGxldCB0aW1lSW5kZXggPSAtMTtcbiAgICBmb3IgKGxldCBpID0gMDsgaSA8IGZyYW1lc1Jvd3MubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgaWYgKGZyYW1lc1Jvd3NbaV0udG9Mb3dlckNhc2UoKS5pbmNsdWRlcyhcInRpbWVcIikpIHtcbiAgICAgICAgICAgIHRpbWVJbmRleCA9IGk7XG4gICAgICAgICAgICBicmVhaztcbiAgICAgICAgfVxuICAgIH1cbiAgICBpZiAodGltZUluZGV4IDwgMCkge1xuICAgICAgICByZXR1cm4geyBmcmFtZXNMZW5ndGg6IDAsIGZyYW1lRHVyYXRpb246IDAsIGZyYW1lczogW10gfTtcbiAgICB9XG4gICAgY29uc3QgZnJhbWVzTGVuZ3RoID0gcGFyc2VJbnQoZnJhbWVzUm93c1t0aW1lSW5kZXggLSAxXS5zcGxpdChcIiBcIikubWFwKChpdGVtKSA9PiBpdGVtLnRyaW0oKSkuZmlsdGVyKChpdGVtKSA9PiAhIWl0ZW0pLmZpbHRlcigoaXRlbSkgPT4gIWlzTmFOKGl0ZW0pKVswXSk7XG4gICAgY29uc3QgZnJhbWVEdXJhdGlvbiA9IHBhcnNlRmxvYXQoZnJhbWVzUm93c1t0aW1lSW5kZXhdLnNwbGl0KFwiIFwiKS5tYXAoKGl0ZW0pID0+IGl0ZW0udHJpbSgpKS5maWx0ZXIoKGl0ZW0pID0+ICEhaXRlbSkuZmlsdGVyKChpdGVtKSA9PiAhaXNOYU4oaXRlbSkpWzBdKTtcbiAgICB3aGlsZSAoIWZyYW1lc1Jvd3NbMF0udG9Mb3dlckNhc2UoKS5pbmNsdWRlcyhcInRpbWVcIikpIHtcbiAgICAgICAgZnJhbWVzUm93cy5zaGlmdCgpO1xuICAgIH1cbiAgICBmcmFtZXNSb3dzLnNoaWZ0KCk7XG4gICAgY29uc3QgZnJhbWVzID0gcGFyc2VGcmFtZXMoZnJhbWVzUm93cyk7XG4gICAgcmV0dXJuIHsgZnJhbWVzTGVuZ3RoLCBmcmFtZUR1cmF0aW9uLCBmcmFtZXMgfTtcbn1cbmZ1bmN0aW9uIGRpc3RyaWJ1dGVTaW5nbGVGcmFtZShoaWVyYXJjaHksIGZyYW1lKSB7XG4gICAgdmFyIF9hLCBfYjtcbiAgICAoX2EgPSBoaWVyYXJjaHkuY2hpbGRyZW4pID09PSBudWxsIHx8IF9hID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfYS50b1JldmVyc2VkKCkuZm9yRWFjaCgoY2hpbGQpID0+IGRpc3RyaWJ1dGVTaW5nbGVGcmFtZShjaGlsZCwgZnJhbWUpKTtcbiAgICBjb25zdCBwb3NpdGlvbiA9IHsgeDogMCwgeTogMCwgejogMCB9O1xuICAgIGNvbnN0IHJvdGF0aW9uID0geyB4OiAwLCB5OiAwLCB6OiAwIH07XG4gICAgKF9iID0gaGllcmFyY2h5LmNoYW5uZWxzKSA9PT0gbnVsbCB8fCBfYiA9PT0gdm9pZCAwID8gdm9pZCAwIDogX2IudG9SZXZlcnNlZCgpLmZvckVhY2goKGl0ZW0pID0+IHtcbiAgICAgICAgY29uc3QgdmFsdWUgPSBmcmFtZS5wb3AoKTtcbiAgICAgICAgKDAsIHV0aWxzXzEuZGlzdHJpYnV0ZVZhbHVlKShwb3NpdGlvbiwgcm90YXRpb24sIGl0ZW0sIHZhbHVlKTtcbiAgICB9KTtcbiAgICBpZiAoIWhpZXJhcmNoeS5idmhGcmFtZXMpIHtcbiAgICAgICAgaGllcmFyY2h5LmJ2aEZyYW1lcyA9IFtdO1xuICAgIH1cbiAgICBoaWVyYXJjaHkuYnZoRnJhbWVzLnB1c2goeyBwb3NpdGlvbiwgcm90YXRpb24gfSk7XG59XG5mdW5jdGlvbiBwYXJzZUJWSCh0ZXh0KSB7XG4gICAgY29uc3QgcGFydHMgPSB0ZXh0LnJlcGxhY2VBbGwoXCJcXHJcIiwgXCJcIikucmVwbGFjZUFsbChcIlxcdFwiLCBcIiBcIikuc3BsaXQoXCJNT1RJT05cIik7XG4gICAgbGV0IHJlc3VsdCA9IHBhcnRzWzBdLnNwbGl0KFwiSElFUkFSQ0hZXCIpWzFdO1xuICAgIFwiYWJjZGVmZ2hpamtsbW5vcHFyc3R1dnd4eXpBQkNERUZHSElKS0xNTk9QUVJTVFVWV1hZWjAxMjM0NTY3ODkgXCIuc3BsaXQoXCJcIikuZm9yRWFjaChpdGVtID0+IHtcbiAgICAgICAgcmVzdWx0ID0gcmVzdWx0LnJlcGxhY2VBbGwoaXRlbSArIFwiXFxuXCIsIGl0ZW0gKyBcIlxcXCIsXFxuXCIpO1xuICAgIH0pO1xuICAgIHJlc3VsdCA9IHJlc3VsdC5yZXBsYWNlQWxsKFwiSk9JTlRcIiwgXCJcXFwiSk9JTlRcXFwiOlxcXCJcIik7XG4gICAgcmVzdWx0ID0gcmVzdWx0LnJlcGxhY2VBbGwoXCJPRkZTRVRcIiwgXCJcXFwiT0ZGU0VUXFxcIjpcXFwiXCIpO1xuICAgIHJlc3VsdCA9IHJlc3VsdC5yZXBsYWNlQWxsKFwiQ0hBTk5FTFNcIiwgXCJcXFwiQ0hBTk5FTFNcXFwiOlxcXCJcIik7XG4gICAgcmVzdWx0ID0gcmVzdWx0LnJlcGxhY2VBbGwoXCJST09UXCIsIFwiXFxcIlJPT1RcXFwiOlxcXCJcIik7XG4gICAgcmVzdWx0ID0gcmVzdWx0LnJlcGxhY2VBbGwoXCJFbmQgU2l0ZVwiLCBcIlxcXCJFbmQgU2l0ZVxcXCI6XFxcIlwiKTtcbiAgICByZXN1bHQgPSByZXN1bHQucmVwbGFjZUFsbChcIntcIiwgXCJcXFwiY29udGVudFxcXCI6IHtcIik7XG4gICAgcmVzdWx0ID0gcmVzdWx0LnNwbGl0KFwifVwiKS5tYXAoaXRlbSA9PiB7XG4gICAgICAgIGlmIChpdGVtLnRyaW0oKS5lbmRzV2l0aChcIixcIikpIHtcbiAgICAgICAgICAgIHJldHVybiBpdGVtLnRyaW0oKSArIFwiXFxcImR1bW15XFxcIjoge31cIjtcbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gaXRlbTtcbiAgICB9KS5qb2luKFwifVwiKTtcbiAgICByZXN1bHQgPSByZXN1bHQuc3BsaXQoXCJ9XCIpLm1hcChpdGVtID0+IHtcbiAgICAgICAgaWYgKGl0ZW0udHJpbSgpLnN0YXJ0c1dpdGgoXCJcXFwiSk9JTlRcXFwiXCIpKSB7XG4gICAgICAgICAgICByZXR1cm4gXCIsXCIgKyBpdGVtLnRyaW0oKTtcbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gaXRlbTtcbiAgICB9KS5qb2luKFwifVwiKTtcbiAgICBsZXQgY291bnQgPSAwO1xuICAgIHJlc3VsdCA9IGVudW1lcmF0ZShyZXN1bHQsIFwiSk9JTlRcIiwgXCJqbnRcIik7XG4gICAgcmVzdWx0ID0gZW51bWVyYXRlKHJlc3VsdCwgXCJjb250ZW50XCIsIFwiY250XCIpO1xuICAgIGNvbnN0IGhpZXJhcmNoeSA9IGNsZWFudXAoSlNPTi5wYXJzZShcIntcIiArIHJlc3VsdCArIFwifVwiKSkuY2hpbGRyZW5bMF07XG4gICAgY29uc3QgYW5pbWF0aW9uID0gcGFyc2VGcmFtZXNQYXJ0KHBhcnRzWzFdKTtcbiAgICBoaWVyYXJjaHkuYnZoVGltZXMgPSBbXTtcbiAgICBmb3IgKGxldCBpID0gMDsgaSA8IGFuaW1hdGlvbi5mcmFtZXNMZW5ndGg7IGkrKykge1xuICAgICAgICBoaWVyYXJjaHkuYnZoVGltZXMucHVzaChhbmltYXRpb24uZnJhbWVEdXJhdGlvbiAqIGkpO1xuICAgIH1cbiAgICBhbmltYXRpb24uZnJhbWVzLmZvckVhY2goaXRlbSA9PiBkaXN0cmlidXRlU2luZ2xlRnJhbWUoaGllcmFyY2h5LCBpdGVtKSk7XG4gICAgcmV0dXJuIGhpZXJhcmNoeTtcbn1cbiIsIlwidXNlIHN0cmljdFwiO1xuT2JqZWN0LmRlZmluZVByb3BlcnR5KGV4cG9ydHMsIFwiX19lc01vZHVsZVwiLCB7IHZhbHVlOiB0cnVlIH0pO1xuZXhwb3J0cy5SQURfVE9fREVHID0gdm9pZCAwO1xuZXhwb3J0cy50b1F1YXRlcm5pb24gPSB0b1F1YXRlcm5pb247XG5leHBvcnRzLnRvRXVsZXJzID0gdG9FdWxlcnM7XG5leHBvcnRzLnF1YXRlcm5pb25Ub0V1bGVycyA9IHF1YXRlcm5pb25Ub0V1bGVycztcbmV4cG9ydHMuYXBwZW5kID0gYXBwZW5kO1xuZXhwb3J0cy5mbG9hdFRvU3RyaW5nID0gZmxvYXRUb1N0cmluZztcbmV4cG9ydHMuZ2V0VW5pZm9ybVRpbWVzID0gZ2V0VW5pZm9ybVRpbWVzO1xuZXhwb3J0cy5jbGlwVGltZXNUb0Nsb3Nlc3RCVkhUaW1lID0gY2xpcFRpbWVzVG9DbG9zZXN0QlZIVGltZTtcbmV4cG9ydHMuZGlzdHJpYnV0ZVZhbHVlID0gZGlzdHJpYnV0ZVZhbHVlO1xuZXhwb3J0cy5sZXJwVmVjdG9yID0gbGVycFZlY3RvcjtcbmV4cG9ydHMubGVycFF1YXRlcm5pb24gPSBsZXJwUXVhdGVybmlvbjtcbmV4cG9ydHMubGVycFZhbHVlcyA9IGxlcnBWYWx1ZXM7XG5jb25zdCBxdWF0ZXJuaW9uXzEgPSByZXF1aXJlKFwicXVhdGVybmlvblwiKTtcbmV4cG9ydHMuUkFEX1RPX0RFRyA9IDE4MC4wIC8gTWF0aC5QSTtcbmZ1bmN0aW9uIHRvUXVhdGVybmlvbih2KSB7XG4gICAgY29uc3Qgd1NxciA9IDEuMCAtICh2LnggKiB2LnggKyB2LnkgKiB2LnkgKyB2LnogKiB2LnopO1xuICAgIHJldHVybiBuZXcgcXVhdGVybmlvbl8xLlF1YXRlcm5pb24oeyB4OiB2LngsIHk6IHYueSwgejogdi56LCB3OiB3U3FyID4gMCA/IE1hdGguc3FydCh3U3FyKSA6IDAgfSk7XG59XG5mdW5jdGlvbiB0b0V1bGVycyh0cnVuY2F0ZWRRdWFudGVyaW9uKSB7XG4gICAgcmV0dXJuIHRvUXVhdGVybmlvbih7XG4gICAgICAgIHg6IHRydW5jYXRlZFF1YW50ZXJpb24ueixcbiAgICAgICAgeTogdHJ1bmNhdGVkUXVhbnRlcmlvbi54LFxuICAgICAgICB6OiB0cnVuY2F0ZWRRdWFudGVyaW9uLnlcbiAgICB9KS50b0V1bGVyKCkubWFwKGl0ZW0xID0+IGl0ZW0xICogMTgwIC8gTWF0aC5QSSk7XG59XG5mdW5jdGlvbiBxdWF0ZXJuaW9uVG9FdWxlcnMocXVhdGVybmlvbikge1xuICAgIGNvbnN0IGV1bGVycyA9IG5ldyBxdWF0ZXJuaW9uXzEuUXVhdGVybmlvbih7XG4gICAgICAgIHg6IHF1YXRlcm5pb24ueixcbiAgICAgICAgeTogcXVhdGVybmlvbi54LFxuICAgICAgICB6OiBxdWF0ZXJuaW9uLnksXG4gICAgICAgIHc6IHF1YXRlcm5pb24ud1xuICAgIH0pLnRvRXVsZXIoKS5tYXAoaXRlbTEgPT4gaXRlbTEgKiAxODAgLyBNYXRoLlBJKTtcbiAgICByZXR1cm4geyB4OiBldWxlcnNbMF0sIHk6IGV1bGVyc1sxXSwgejogZXVsZXJzWzJdIH07XG59XG5mdW5jdGlvbiBhcHBlbmQodGV4dCwgdGltZXMpIHtcbiAgICBsZXQgcmVzdWx0ID0gXCJcIjtcbiAgICBmb3IgKGxldCBpID0gMDsgaSA8IHRpbWVzOyBpKyspIHtcbiAgICAgICAgcmVzdWx0ICs9IHRleHQ7XG4gICAgfVxuICAgIHJldHVybiByZXN1bHQ7XG59XG5mdW5jdGlvbiBmbG9hdFRvU3RyaW5nKHZhbHVlLCBmcmFjdGlvbikge1xuICAgIHJldHVybiB2YWx1ZS50b0xvY2FsZVN0cmluZyhcInVuLVVTXCIsIHsgbWluaW11bUZyYWN0aW9uRGlnaXRzOiBmcmFjdGlvbiB9KS5yZXBsYWNlQWxsKFwiLFwiLCBcIlwiKTtcbn1cbmZ1bmN0aW9uIGxlcnBWYWx1ZSh4MSwgeDIsIHQpIHtcbiAgICBpZiAodCA+IDEpIHtcbiAgICAgICAgcmV0dXJuIHgyO1xuICAgIH1cbiAgICBpZiAodCA8IDApIHtcbiAgICAgICAgcmV0dXJuIHgxO1xuICAgIH1cbiAgICByZXR1cm4geDEgKyAoeDIgLSB4MSkgKiB0O1xufVxuZnVuY3Rpb24gb3B0aW1pemVkQW1pbWF0aW9uTGVuZ3RoKGR1cmF0aW9uLCBvcmlnaW5hbEZyYW1lTGVuZ3RoKSB7XG4gICAgaWYgKG9yaWdpbmFsRnJhbWVMZW5ndGggPiBkdXJhdGlvbikge1xuICAgICAgICByZXR1cm4gMjtcbiAgICB9XG4gICAgY29uc3QgaG91ckxlbmd0aCA9IDM2MDAgLyBvcmlnaW5hbEZyYW1lTGVuZ3RoO1xuICAgIGxldCBiZXN0TGVuZ3RoID0gMDtcbiAgICBmb3IgKGxldCBpID0gMTsgaSA8IGhvdXJMZW5ndGg7IGkrKykge1xuICAgICAgICBjb25zdCBvcHRpbWl6ZWRGcmFtZUxlbmd0aCA9IGR1cmF0aW9uIC8gaTtcbiAgICAgICAgY29uc3QgZXJyb3IgPSBNYXRoLmFicyhvcmlnaW5hbEZyYW1lTGVuZ3RoIC0gb3B0aW1pemVkRnJhbWVMZW5ndGgpO1xuICAgICAgICBjb25zdCBwcmV2RXJyb3IgPSBNYXRoLmFicyhvcmlnaW5hbEZyYW1lTGVuZ3RoIC0gYmVzdExlbmd0aCk7XG4gICAgICAgIGlmIChlcnJvciA8IHByZXZFcnJvcikge1xuICAgICAgICAgICAgYmVzdExlbmd0aCA9IG9wdGltaXplZEZyYW1lTGVuZ3RoO1xuICAgICAgICB9XG4gICAgfVxuICAgIHJldHVybiBiZXN0TGVuZ3RoO1xufVxuZnVuY3Rpb24gZ2V0VW5pZm9ybVRpbWVzKGR1cmF0aW9uLCBzaW5nbGVGcmFtZUR1cmF0aW9uKSB7XG4gICAgY29uc3QgdGltZXMgPSBbXTtcbiAgICBjb25zdCBvcHRpbWl6ZWRGcmFtZUR1cmF0aW9uID0gb3B0aW1pemVkQW1pbWF0aW9uTGVuZ3RoKGR1cmF0aW9uLCBzaW5nbGVGcmFtZUR1cmF0aW9uKTtcbiAgICBjb25zdCBsZW5ndGggPSBkdXJhdGlvbiAvIG9wdGltaXplZEZyYW1lRHVyYXRpb247XG4gICAgZm9yIChsZXQgaSA9IDA7IGkgPCBsZW5ndGg7IGkrKykge1xuICAgICAgICB0aW1lcy5wdXNoKGkgKiBvcHRpbWl6ZWRGcmFtZUR1cmF0aW9uKTtcbiAgICB9XG4gICAgdGltZXNbbGVuZ3RoIC0gMV0gPSBkdXJhdGlvbjtcbiAgICByZXR1cm4gdGltZXM7XG59XG5mdW5jdGlvbiBjbG9zZXN0KGxlZnQsIHJpZ2h0LCB2YWx1ZSkge1xuICAgIGlmIChNYXRoLmFicyhsZWZ0IC0gdmFsdWUpIDwgTWF0aC5hYnMocmlnaHQgLSB2YWx1ZSkpIHtcbiAgICAgICAgcmV0dXJuIGxlZnQ7XG4gICAgfVxuICAgIHJldHVybiByaWdodDtcbn1cbmZ1bmN0aW9uIGNsaXBUaW1lc1RvQ2xvc2VzdEJWSFRpbWUoYW5pbVRpbWVzLCBidmhUaW1lcykge1xuICAgIGNvbnN0IGZpeGVkVGltZXMgPSBhbmltVGltZXMubWFwKChpdGVtKSA9PiBpdGVtKTtcbiAgICBmb3IgKGxldCBpID0gMTsgaSA8IGJ2aFRpbWVzLmxlbmd0aDsgaSsrKSB7XG4gICAgICAgIGZvciAobGV0IGogPSAwOyBqIDwgYW5pbVRpbWVzLmxlbmd0aDsgaisrKSB7XG4gICAgICAgICAgICBjb25zdCBhbmltVGltZSA9IGFuaW1UaW1lc1tqXTtcbiAgICAgICAgICAgIGNvbnN0IGJ2aFRpbWVMZWZ0ID0gYnZoVGltZXNbaSAtIDFdO1xuICAgICAgICAgICAgY29uc3QgYnZoVGltZVJpZ2h0ID0gYnZoVGltZXNbaV07XG4gICAgICAgICAgICBpZiAoKGJ2aFRpbWVMZWZ0IDw9IGFuaW1UaW1lKSAmJiAoYW5pbVRpbWUgPD0gYnZoVGltZVJpZ2h0KSkge1xuICAgICAgICAgICAgICAgIGZpeGVkVGltZXNbal0gPSBjbG9zZXN0KGJ2aFRpbWVMZWZ0LCBidmhUaW1lUmlnaHQsIGFuaW1UaW1lKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgIH1cbiAgICByZXR1cm4gZml4ZWRUaW1lcztcbn1cbmZ1bmN0aW9uIGdldEZhY3RvcnMoYnZoVGltZXMsIGFuaW1UaW1lcykge1xuICAgIHJldHVybiBidmhUaW1lcy5tYXAoKGl0ZW0pID0+IHtcbiAgICAgICAgaWYgKGl0ZW0gPD0gYW5pbVRpbWVzWzBdKSB7XG4gICAgICAgICAgICByZXR1cm4ge1xuICAgICAgICAgICAgICAgIGxlZnRBbmltSW5kZXg6IDAsXG4gICAgICAgICAgICAgICAgcmlnaHRBbmltSW5kZXg6IDEsXG4gICAgICAgICAgICAgICAgZmFjdG9yOiAwXG4gICAgICAgICAgICB9O1xuICAgICAgICB9XG4gICAgICAgIGlmIChpdGVtID49IGFuaW1UaW1lc1thbmltVGltZXMubGVuZ3RoIC0gMV0pIHtcbiAgICAgICAgICAgIHJldHVybiB7XG4gICAgICAgICAgICAgICAgbGVmdEFuaW1JbmRleDogYW5pbVRpbWVzLmxlbmd0aCAtIDIsXG4gICAgICAgICAgICAgICAgcmlnaHRBbmltSW5kZXg6IGFuaW1UaW1lcy5sZW5ndGggLSAxLFxuICAgICAgICAgICAgICAgIGZhY3RvcjogMVxuICAgICAgICAgICAgfTtcbiAgICAgICAgfVxuICAgICAgICBmb3IgKGxldCBpID0gMTsgaSA8IGFuaW1UaW1lcy5sZW5ndGg7IGkrKykge1xuICAgICAgICAgICAgY29uc3QgbGVmdFRpbWUgPSBhbmltVGltZXNbaSAtIDFdO1xuICAgICAgICAgICAgY29uc3QgcmlnaHRUaW1lID0gYW5pbVRpbWVzW2ldO1xuICAgICAgICAgICAgaWYgKChsZWZ0VGltZSA8PSBpdGVtKSAmJiAoaXRlbSA8IHJpZ2h0VGltZSkpIHtcbiAgICAgICAgICAgICAgICBjb25zdCByYW5nZVNpemUgPSByaWdodFRpbWUgLSBsZWZ0VGltZTtcbiAgICAgICAgICAgICAgICByZXR1cm4ge1xuICAgICAgICAgICAgICAgICAgICBsZWZ0QW5pbUluZGV4OiBpIC0gMSxcbiAgICAgICAgICAgICAgICAgICAgcmlnaHRBbmltSW5kZXg6IGksXG4gICAgICAgICAgICAgICAgICAgIGZhY3RvcjogKGl0ZW0gLSBsZWZ0VGltZSkgLyByYW5nZVNpemVcbiAgICAgICAgICAgICAgICB9O1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgfSk7XG59XG5mdW5jdGlvbiBkaXN0cmlidXRlVmFsdWUocG9zaXRpb24sIHJvdGF0aW9uLCBjaGFubmVsLCB2YWx1ZSkge1xuICAgIGNvbnN0IGtleSA9IGNoYW5uZWwudG9Mb3dlckNhc2UoKVswXTtcbiAgICBjb25zdCByZWNpcGllbnQgPSBjaGFubmVsLmluY2x1ZGVzKFwicG9zXCIpID8gcG9zaXRpb24gOiByb3RhdGlvbjtcbiAgICByZWNpcGllbnRba2V5XSA9IHZhbHVlO1xufVxuZnVuY3Rpb24gbGVycFZlY3RvcihsZWZ0VmFsdWUsIHJpZ2h0VmFsdWUsIGZhY3Rvcikge1xuICAgIHJldHVybiB7XG4gICAgICAgIHg6IGxlcnBWYWx1ZShsZWZ0VmFsdWUueCwgcmlnaHRWYWx1ZS54LCBmYWN0b3IpLFxuICAgICAgICB5OiBsZXJwVmFsdWUobGVmdFZhbHVlLnksIHJpZ2h0VmFsdWUueSwgZmFjdG9yKSxcbiAgICAgICAgejogbGVycFZhbHVlKGxlZnRWYWx1ZS56LCByaWdodFZhbHVlLnosIGZhY3RvciksXG4gICAgfTtcbn1cbmZ1bmN0aW9uIGxlcnBRdWF0ZXJuaW9uKGxlZnRWYWx1ZSwgcmlnaHRWYWx1ZSwgZmFjdG9yKSB7XG4gICAgcmV0dXJuIGxlZnRWYWx1ZS5zbGVycChyaWdodFZhbHVlKShmYWN0b3IpO1xufVxuZnVuY3Rpb24gbGVycFZhbHVlcyh2YWx1ZXMsIGFuaW1UaW1lcywgdW5pZm9ybVRpbWVzLCBkZWZhdWx0VmFsdWUsIGxlcnBGdW5jdGlvbikge1xuICAgIGlmICghKHZhbHVlcyA9PT0gbnVsbCB8fCB2YWx1ZXMgPT09IHZvaWQgMCA/IHZvaWQgMCA6IHZhbHVlcy5sZW5ndGgpKSB7XG4gICAgICAgIHJldHVybiB1bmlmb3JtVGltZXMubWFwKGl0ZW0gPT4gZGVmYXVsdFZhbHVlKTtcbiAgICB9XG4gICAgaWYgKHZhbHVlcy5sZW5ndGggPT0gMSkge1xuICAgICAgICByZXR1cm4gdW5pZm9ybVRpbWVzLm1hcChpdGVtID0+IHZhbHVlc1swXSk7XG4gICAgfVxuICAgIGNvbnN0IGZhY3RvcnMgPSBnZXRGYWN0b3JzKHVuaWZvcm1UaW1lcywgYW5pbVRpbWVzKTtcbiAgICByZXR1cm4gZmFjdG9ycy5tYXAoaXRlbSA9PiB7XG4gICAgICAgIGNvbnN0IGxlZnRGcmFtZSA9IHZhbHVlc1tpdGVtLmxlZnRBbmltSW5kZXhdO1xuICAgICAgICBjb25zdCByaWdodEZyYW1lID0gdmFsdWVzW2l0ZW0ucmlnaHRBbmltSW5kZXhdO1xuICAgICAgICByZXR1cm4gbGVycEZ1bmN0aW9uKGxlZnRGcmFtZSwgcmlnaHRGcmFtZSwgaXRlbS5mYWN0b3IpO1xuICAgIH0pO1xufVxuIiwiLy8gVGhlIG1vZHVsZSBjYWNoZVxudmFyIF9fd2VicGFja19tb2R1bGVfY2FjaGVfXyA9IHt9O1xuXG4vLyBUaGUgcmVxdWlyZSBmdW5jdGlvblxuZnVuY3Rpb24gX193ZWJwYWNrX3JlcXVpcmVfXyhtb2R1bGVJZCkge1xuXHQvLyBDaGVjayBpZiBtb2R1bGUgaXMgaW4gY2FjaGVcblx0dmFyIGNhY2hlZE1vZHVsZSA9IF9fd2VicGFja19tb2R1bGVfY2FjaGVfX1ttb2R1bGVJZF07XG5cdGlmIChjYWNoZWRNb2R1bGUgIT09IHVuZGVmaW5lZCkge1xuXHRcdHJldHVybiBjYWNoZWRNb2R1bGUuZXhwb3J0cztcblx0fVxuXHQvLyBDcmVhdGUgYSBuZXcgbW9kdWxlIChhbmQgcHV0IGl0IGludG8gdGhlIGNhY2hlKVxuXHR2YXIgbW9kdWxlID0gX193ZWJwYWNrX21vZHVsZV9jYWNoZV9fW21vZHVsZUlkXSA9IHtcblx0XHQvLyBubyBtb2R1bGUuaWQgbmVlZGVkXG5cdFx0Ly8gbm8gbW9kdWxlLmxvYWRlZCBuZWVkZWRcblx0XHRleHBvcnRzOiB7fVxuXHR9O1xuXG5cdC8vIEV4ZWN1dGUgdGhlIG1vZHVsZSBmdW5jdGlvblxuXHRfX3dlYnBhY2tfbW9kdWxlc19fW21vZHVsZUlkXShtb2R1bGUsIG1vZHVsZS5leHBvcnRzLCBfX3dlYnBhY2tfcmVxdWlyZV9fKTtcblxuXHQvLyBSZXR1cm4gdGhlIGV4cG9ydHMgb2YgdGhlIG1vZHVsZVxuXHRyZXR1cm4gbW9kdWxlLmV4cG9ydHM7XG59XG5cbiIsIlwidXNlIHN0cmljdFwiO1xuT2JqZWN0LmRlZmluZVByb3BlcnR5KGV4cG9ydHMsIFwiX19lc01vZHVsZVwiLCB7IHZhbHVlOiB0cnVlIH0pO1xuZXhwb3J0cy5mZW1hbGVPZmZzZXRzID0gZXhwb3J0cy5tYWxlT2Zmc2V0cyA9IGV4cG9ydHMudmlzaXROb2RlID0gZXhwb3J0cy5zZXJpYWxpemVCVkggPSBleHBvcnRzLnRvQlZIID0gZXhwb3J0cy5wYXJzZUJWSCA9IGV4cG9ydHMucGFyc2VBbmltID0gdm9pZCAwO1xuY29uc3QgcGFyc2VfMSA9IHJlcXVpcmUoXCIuL3BhcnNlXCIpO1xuT2JqZWN0LmRlZmluZVByb3BlcnR5KGV4cG9ydHMsIFwicGFyc2VBbmltXCIsIHsgZW51bWVyYWJsZTogdHJ1ZSwgZ2V0OiBmdW5jdGlvbiAoKSB7IHJldHVybiBwYXJzZV8xLnBhcnNlQW5pbTsgfSB9KTtcbk9iamVjdC5kZWZpbmVQcm9wZXJ0eShleHBvcnRzLCBcInBhcnNlQlZIXCIsIHsgZW51bWVyYWJsZTogdHJ1ZSwgZ2V0OiBmdW5jdGlvbiAoKSB7IHJldHVybiBwYXJzZV8xLnBhcnNlQlZIOyB9IH0pO1xuY29uc3QgY29udmVydF8xID0gcmVxdWlyZShcIi4vY29udmVydFwiKTtcbk9iamVjdC5kZWZpbmVQcm9wZXJ0eShleHBvcnRzLCBcInRvQlZIXCIsIHsgZW51bWVyYWJsZTogdHJ1ZSwgZ2V0OiBmdW5jdGlvbiAoKSB7IHJldHVybiBjb252ZXJ0XzEudG9CVkg7IH0gfSk7XG5PYmplY3QuZGVmaW5lUHJvcGVydHkoZXhwb3J0cywgXCJzZXJpYWxpemVCVkhcIiwgeyBlbnVtZXJhYmxlOiB0cnVlLCBnZXQ6IGZ1bmN0aW9uICgpIHsgcmV0dXJuIGNvbnZlcnRfMS5zZXJpYWxpemVCVkg7IH0gfSk7XG5PYmplY3QuZGVmaW5lUHJvcGVydHkoZXhwb3J0cywgXCJ2aXNpdE5vZGVcIiwgeyBlbnVtZXJhYmxlOiB0cnVlLCBnZXQ6IGZ1bmN0aW9uICgpIHsgcmV0dXJuIGNvbnZlcnRfMS52aXNpdE5vZGU7IH0gfSk7XG5jb25zdCBvZmZzZXRzXzEgPSByZXF1aXJlKFwiLi9vZmZzZXRzXCIpO1xuT2JqZWN0LmRlZmluZVByb3BlcnR5KGV4cG9ydHMsIFwibWFsZU9mZnNldHNcIiwgeyBlbnVtZXJhYmxlOiB0cnVlLCBnZXQ6IGZ1bmN0aW9uICgpIHsgcmV0dXJuIG9mZnNldHNfMS5tYWxlT2Zmc2V0czsgfSB9KTtcbk9iamVjdC5kZWZpbmVQcm9wZXJ0eShleHBvcnRzLCBcImZlbWFsZU9mZnNldHNcIiwgeyBlbnVtZXJhYmxlOiB0cnVlLCBnZXQ6IGZ1bmN0aW9uICgpIHsgcmV0dXJuIG9mZnNldHNfMS5mZW1hbGVPZmZzZXRzOyB9IH0pO1xuIl0sIm5hbWVzIjpbXSwic291cmNlUm9vdCI6IiJ9