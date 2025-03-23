(function webpackUniversalModuleDefinition(root, factory) {
	if(typeof exports === 'object' && typeof module === 'object')
		module.exports = factory();
	else if(typeof define === 'function' && define.amd)
		define([], factory);
	else if(typeof exports === 'object')
		exports["AnimToBvh"] = factory();
	else
		root["AnimToBvh"] = factory();
})(self, () => {
return /******/ (() => { // webpackBootstrap
/******/ 	"use strict";
/******/ 	var __webpack_modules__ = ({

/***/ "./node_modules/quaternion/dist/quaternion.mjs":
/*!*****************************************************!*\
  !*** ./node_modules/quaternion/dist/quaternion.mjs ***!
  \*****************************************************/
/***/ ((__unused_webpack___webpack_module__, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Quaternion: () => (/* binding */ Quaternion),
/* harmony export */   "default": () => (/* binding */ Quaternion)
/* harmony export */ });


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



/***/ }),

/***/ "./src/aliases.ts":
/*!************************!*\
  !*** ./src/aliases.ts ***!
  \************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   aliases: () => (/* binding */ aliases)
/* harmony export */ });
const aliases = {
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
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   serializeBVH: () => (/* binding */ serializeBVH),
/* harmony export */   toBVH: () => (/* binding */ toBVH),
/* harmony export */   visitNode: () => (/* binding */ visitNode)
/* harmony export */ });
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./utils */ "./src/utils.ts");
/* harmony import */ var _hierarchy__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./hierarchy */ "./src/hierarchy.ts");
/* harmony import */ var _aliases__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./aliases */ "./src/aliases.ts");



function offsetToString(offset, digits) {
    return "OFFSET " + (0,_utils__WEBPACK_IMPORTED_MODULE_0__.floatToString)(offset.x, digits) + " " + (0,_utils__WEBPACK_IMPORTED_MODULE_0__.floatToString)(offset.y, digits) + " " + (0,_utils__WEBPACK_IMPORTED_MODULE_0__.floatToString)(offset.z, digits);
}
function appendNode(joint, tabs) {
    let result = "";
    const boneType = (joint.bvhName === "hip") ? "ROOT" : "JOINT";
    const channels = (joint.bvhName === "hip") ? "CHANNELS 6 Xposition Yposition Zposition Xrotation Yrotation Zrotation" : "CHANNELS 3 Xrotation Yrotation Zrotation";
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
    const bvhNames = names.map(item => _aliases__WEBPACK_IMPORTED_MODULE_2__.aliases[item] || item);
    return collectNodes(_hierarchy__WEBPACK_IMPORTED_MODULE_1__.hierarchy, bvhNames);
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
function fillKeyFrames(data, bvhNode, fps) {
    const animJoints = data.joints;
    const length = extractFramesLength(animJoints);
    const animTimes = extractTimes(animJoints);
    const bvhTimes = (0,_utils__WEBPACK_IMPORTED_MODULE_0__.getUniformTimes)(data.duration, 1 / fps);
    const fixedAnimTimes = (0,_utils__WEBPACK_IMPORTED_MODULE_0__.clipTimesToClosestBVHTime)(animTimes, bvhTimes);
    visitNode(bvhNode, (node) => {
        const joint = animJoints.find((item) => _aliases__WEBPACK_IMPORTED_MODULE_2__.aliases[item.joint_name] === node.bvhName);
        node.offset = { x: 0, y: 0, z: 0 };
        if (node.bvhName != "end") {
            fillChannels(node, joint);
            node.animFrames = [];
            for (let i = 0; i < length; i++) {
                node.animFrames.push({
                    position: (joint === null || joint === void 0 ? void 0 : joint.position_keys[i]) || { x: 0, y: 0, z: 0 },
                    rotation: (joint === null || joint === void 0 ? void 0 : joint.rotation_keys[i]) || { x: 0, y: 0, z: 0 },
                    time: animTimes[i]
                });
            }
            const positions = (0,_utils__WEBPACK_IMPORTED_MODULE_0__.lerpValues)(node.animFrames.map((item) => item.position), fixedAnimTimes, bvhTimes, _utils__WEBPACK_IMPORTED_MODULE_0__.lerpVector);
            const rotations = (0,_utils__WEBPACK_IMPORTED_MODULE_0__.lerpValues)(node.animFrames.map((item) => (0,_utils__WEBPACK_IMPORTED_MODULE_0__.toQuaternion)(item.rotation)), fixedAnimTimes, bvhTimes, _utils__WEBPACK_IMPORTED_MODULE_0__.lerpQuaternion).map(item => (0,_utils__WEBPACK_IMPORTED_MODULE_0__.quaternionToEulers)(item));
            node.bvhFrames = [];
            positions.forEach((item, i) => node.bvhFrames.push({
                position: item,
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
    return values.map(item => (0,_utils__WEBPACK_IMPORTED_MODULE_0__.floatToString)(item, 4)).join(" ") + " \n";
}
function serializeBVH(bvhNode) {
    let result = "HIERARCHY\n";
    result += appendNode(bvhNode, "");
    result += "MOTION\n";
    result += "Frames " + bvhNode.bvhTimes.length + "\n";
    result += "Frame Time " + (0,_utils__WEBPACK_IMPORTED_MODULE_0__.floatToString)(bvhNode.bvhTimes[1], 6) + "\n";
    for (let i = 0; i < bvhNode.bvhTimes.length; i++) {
        result += getFrameRow(bvhNode, i);
    }
    return result;
}
function toBVH(data) {
    const bvhNode = subTree(data.joints);
    fillKeyFrames(data, bvhNode, 24);
    return bvhNode;
}


/***/ }),

/***/ "./src/hierarchy.ts":
/*!**************************!*\
  !*** ./src/hierarchy.ts ***!
  \**************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   hierarchy: () => (/* binding */ hierarchy)
/* harmony export */ });
const hierarchy = {
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
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   femaleOffsets: () => (/* binding */ femaleOffsets),
/* harmony export */   maleOffsets: () => (/* binding */ maleOffsets)
/* harmony export */ });
const femaleOffsets = {
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
const maleOffsets = {
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
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   distributeSingleFrame: () => (/* binding */ distributeSingleFrame),
/* harmony export */   parseAnim: () => (/* binding */ parseAnim),
/* harmony export */   parseBVH: () => (/* binding */ parseBVH)
/* harmony export */ });
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./utils */ "./src/utils.ts");

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
        item.euler_keys.push((0,_utils__WEBPACK_IMPORTED_MODULE_0__.toEulers)(rot));
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
        const channels = node.CHANNELS.trim().split(" ").filter((item) => isNaN(parseFloat(item)));
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
        (0,_utils__WEBPACK_IMPORTED_MODULE_0__.distributeValue)(position, rotation, item, value);
    });
    if (!hierarchy.bvhFrames) {
        hierarchy.bvhFrames = [];
    }
    hierarchy.bvhFrames.push({ position, rotation });
}
function parseBVH(text) {
    const parts = text.split("MOTION");
    let result = parts[0].split("HIERARCHY")[1];
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789".split("").forEach(item => {
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
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   RAD_TO_DEG: () => (/* binding */ RAD_TO_DEG),
/* harmony export */   append: () => (/* binding */ append),
/* harmony export */   clipTimesToClosestBVHTime: () => (/* binding */ clipTimesToClosestBVHTime),
/* harmony export */   distributeValue: () => (/* binding */ distributeValue),
/* harmony export */   floatToString: () => (/* binding */ floatToString),
/* harmony export */   getUniformTimes: () => (/* binding */ getUniformTimes),
/* harmony export */   lerpQuaternion: () => (/* binding */ lerpQuaternion),
/* harmony export */   lerpValues: () => (/* binding */ lerpValues),
/* harmony export */   lerpVector: () => (/* binding */ lerpVector),
/* harmony export */   quaternionToEulers: () => (/* binding */ quaternionToEulers),
/* harmony export */   toEulers: () => (/* binding */ toEulers),
/* harmony export */   toQuaternion: () => (/* binding */ toQuaternion)
/* harmony export */ });
/* harmony import */ var quaternion__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! quaternion */ "./node_modules/quaternion/dist/quaternion.mjs");

const RAD_TO_DEG = 180.0 / Math.PI;
function toQuaternion(v) {
    const wSqr = 1.0 - (v.x * v.x + v.y * v.y + v.z * v.z);
    return new quaternion__WEBPACK_IMPORTED_MODULE_0__.Quaternion({ x: v.x, y: v.y, z: v.z, w: wSqr > 0 ? Math.sqrt(wSqr) : 0 });
}
function toEulers(truncatedQuanterion) {
    return toQuaternion({
        x: truncatedQuanterion.z,
        y: truncatedQuanterion.x,
        z: truncatedQuanterion.y
    }).toEuler().map(item1 => item1 * 180 / Math.PI);
}
function quaternionToEulers(quaternion) {
    const eulers = new quaternion__WEBPACK_IMPORTED_MODULE_0__.Quaternion({
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
function lerpValues(values, animTimes, uniformTimes, lerpFunction) {
    return getFactors(uniformTimes, animTimes).map(item => {
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
/******/ 	/* webpack/runtime/define property getters */
/******/ 	(() => {
/******/ 		// define getter functions for harmony exports
/******/ 		__webpack_require__.d = (exports, definition) => {
/******/ 			for(var key in definition) {
/******/ 				if(__webpack_require__.o(definition, key) && !__webpack_require__.o(exports, key)) {
/******/ 					Object.defineProperty(exports, key, { enumerable: true, get: definition[key] });
/******/ 				}
/******/ 			}
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/hasOwnProperty shorthand */
/******/ 	(() => {
/******/ 		__webpack_require__.o = (obj, prop) => (Object.prototype.hasOwnProperty.call(obj, prop))
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/make namespace object */
/******/ 	(() => {
/******/ 		// define __esModule on exports
/******/ 		__webpack_require__.r = (exports) => {
/******/ 			if(typeof Symbol !== 'undefined' && Symbol.toStringTag) {
/******/ 				Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
/******/ 			}
/******/ 			Object.defineProperty(exports, '__esModule', { value: true });
/******/ 		};
/******/ 	})();
/******/ 	
/************************************************************************/
var __webpack_exports__ = {};
// This entry needs to be wrapped in an IIFE because it needs to be isolated against other modules in the chunk.
(() => {
/*!**********************!*\
  !*** ./src/index.ts ***!
  \**********************/
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   femaleOffsets: () => (/* reexport safe */ _offsets__WEBPACK_IMPORTED_MODULE_2__.femaleOffsets),
/* harmony export */   maleOffsets: () => (/* reexport safe */ _offsets__WEBPACK_IMPORTED_MODULE_2__.maleOffsets),
/* harmony export */   parseAnim: () => (/* reexport safe */ _parse__WEBPACK_IMPORTED_MODULE_0__.parseAnim),
/* harmony export */   parseBVH: () => (/* reexport safe */ _parse__WEBPACK_IMPORTED_MODULE_0__.parseBVH),
/* harmony export */   serializeBVH: () => (/* reexport safe */ _convert__WEBPACK_IMPORTED_MODULE_1__.serializeBVH),
/* harmony export */   toBVH: () => (/* reexport safe */ _convert__WEBPACK_IMPORTED_MODULE_1__.toBVH),
/* harmony export */   visitNode: () => (/* reexport safe */ _convert__WEBPACK_IMPORTED_MODULE_1__.visitNode)
/* harmony export */ });
/* harmony import */ var _parse__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./parse */ "./src/parse.ts");
/* harmony import */ var _convert__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./convert */ "./src/convert.ts");
/* harmony import */ var _offsets__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./offsets */ "./src/offsets.ts");





})();

/******/ 	return __webpack_exports__;
/******/ })()
;
});
//# sourceMappingURL=data:application/json;charset=utf-8;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoiYnVuZGxlLmpzIiwibWFwcGluZ3MiOiJBQUFBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLENBQUM7QUFDRCxPOzs7Ozs7Ozs7Ozs7Ozs7QUNWYTs7QUFFYjtBQUNBO0FBQ0E7QUFDQSxXQUFXLFFBQVE7QUFDbkIsV0FBVyxRQUFRO0FBQ25CLFdBQVcsUUFBUTtBQUNuQixXQUFXLFFBQVE7QUFDbkI7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQSxXQUFXLFFBQVE7QUFDbkIsV0FBVyxRQUFRO0FBQ25CLFdBQVcsUUFBUTtBQUNuQixXQUFXLFFBQVE7QUFDbkI7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBLFdBQVcsUUFBUTtBQUNuQixXQUFXLFFBQVE7QUFDbkIsYUFBYTtBQUNiO0FBQ0E7O0FBRUE7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBOztBQUVBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTs7QUFFQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQSxtQkFBbUI7O0FBRW5CO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBLG9CQUFvQixtQkFBbUI7O0FBRXZDO0FBQ0E7O0FBRUE7QUFDQTtBQUNBLFFBQVE7QUFDUjtBQUNBLFFBQVE7QUFDUjtBQUNBLFFBQVE7O0FBRVI7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLFlBQVk7QUFDWjtBQUNBOztBQUVBLFVBQVU7O0FBRVY7QUFDQTtBQUNBOztBQUVBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxJQUFJOztBQUVKOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxNQUFNO0FBQ047QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBOztBQUVBOztBQUVBOztBQUVBO0FBQ0E7QUFDQSxNQUFNO0FBQ047QUFDQTs7QUFFQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLFdBQVcsc0JBQXNCO0FBQ2pDLFdBQVcsU0FBUztBQUNwQixXQUFXLFNBQVM7QUFDcEIsV0FBVyxTQUFTO0FBQ3BCLGFBQWE7QUFDYjtBQUNBOztBQUVBO0FBQ0E7QUFDQSxJQUFJO0FBQ0o7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsYUFBYSxzQkFBc0I7QUFDbkMsYUFBYSxTQUFTO0FBQ3RCLGFBQWEsU0FBUztBQUN0QixhQUFhLFNBQVM7QUFDdEIsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7O0FBRUE7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQSxhQUFhLHNCQUFzQjtBQUNuQyxhQUFhLFNBQVM7QUFDdEIsYUFBYSxTQUFTO0FBQ3RCLGFBQWEsU0FBUztBQUN0QixlQUFlO0FBQ2Y7QUFDQTs7QUFFQTs7QUFFQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7O0FBRUE7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7O0FBRUE7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0EsR0FBRztBQUNIO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7O0FBRUE7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7O0FBRUE7QUFDQTtBQUNBOztBQUVBOztBQUVBO0FBQ0EsR0FBRztBQUNIO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsYUFBYSxzQkFBc0I7QUFDbkMsYUFBYSxTQUFTO0FBQ3RCLGFBQWEsU0FBUztBQUN0QixhQUFhLFNBQVM7QUFDdEIsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7O0FBRUE7O0FBRUE7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQSxhQUFhLFFBQVE7QUFDckIsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQSxhQUFhLHNCQUFzQjtBQUNuQyxhQUFhLFNBQVM7QUFDdEIsYUFBYSxTQUFTO0FBQ3RCLGFBQWEsU0FBUztBQUN0QixlQUFlO0FBQ2Y7QUFDQTs7QUFFQTs7QUFFQTs7QUFFQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQTtBQUNBLGVBQWU7QUFDZjtBQUNBOztBQUVBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTs7QUFFQTtBQUNBLGlDQUFpQztBQUNqQzs7QUFFQTs7QUFFQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQSxhQUFhLHNCQUFzQjtBQUNuQyxhQUFhLFNBQVM7QUFDdEIsYUFBYSxTQUFTO0FBQ3RCLGFBQWEsU0FBUztBQUN0QixlQUFlO0FBQ2Y7QUFDQTs7QUFFQTs7QUFFQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTs7QUFFQTtBQUNBLGlDQUFpQztBQUNqQzs7QUFFQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsR0FBRztBQUNIO0FBQ0E7QUFDQTtBQUNBLGVBQWU7QUFDZjtBQUNBOztBQUVBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBOztBQUVBO0FBQ0E7O0FBRUE7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBOztBQUVBLGdEQUFnRDs7QUFFaEQ7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQSxhQUFhLHNCQUFzQjtBQUNuQyxhQUFhLFNBQVM7QUFDdEIsYUFBYSxTQUFTO0FBQ3RCLGFBQWEsU0FBUztBQUN0QixlQUFlO0FBQ2Y7QUFDQTs7QUFFQTs7QUFFQTs7QUFFQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBOztBQUVBOztBQUVBOztBQUVBOztBQUVBLFlBQVk7O0FBRVo7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQSxhQUFhLHNCQUFzQjtBQUNuQyxhQUFhLFNBQVM7QUFDdEIsYUFBYSxTQUFTO0FBQ3RCLGFBQWEsU0FBUztBQUN0QixlQUFlO0FBQ2Y7QUFDQTs7QUFFQTs7QUFFQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsR0FBRztBQUNIO0FBQ0E7QUFDQTtBQUNBLGVBQWU7QUFDZjtBQUNBOztBQUVBO0FBQ0EsR0FBRztBQUNIO0FBQ0E7QUFDQTtBQUNBLGVBQWU7QUFDZjtBQUNBOztBQUVBO0FBQ0EsR0FBRztBQUNIO0FBQ0E7QUFDQTtBQUNBLGVBQWU7QUFDZjtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBOztBQUVBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTs7QUFFQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQSxlQUFlO0FBQ2Y7QUFDQTs7QUFFQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQSxlQUFlO0FBQ2Y7QUFDQTs7QUFFQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQSxlQUFlO0FBQ2Y7QUFDQTs7QUFFQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQSxhQUFhLFVBQVU7QUFDdkIsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQSxhQUFhLFVBQVU7QUFDdkIsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7O0FBRUE7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQSxNQUFNO0FBQ04sa0NBQWtDO0FBQ2xDO0FBQ0E7QUFDQSxHQUFHO0FBQ0g7QUFDQTtBQUNBO0FBQ0EsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7QUFDQSw0QkFBNEI7O0FBRTVCLDBCQUEwQjtBQUMxQjtBQUNBLHFEQUFxRDtBQUNyRDs7QUFFQTtBQUNBLG9DQUFvQztBQUNwQztBQUNBLEdBQUc7QUFDSDtBQUNBO0FBQ0E7QUFDQSxhQUFhLFNBQVM7QUFDdEIsZUFBZTtBQUNmO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQSwrQ0FBK0M7QUFDL0M7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQSxrQ0FBa0M7QUFDbEM7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBLG1DQUFtQztBQUNuQztBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsR0FBRztBQUNIO0FBQ0E7QUFDQTtBQUNBLGVBQWU7QUFDZjtBQUNBOztBQUVBO0FBQ0EsR0FBRztBQUNIO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsYUFBYSxPQUFPO0FBQ3BCLGVBQWU7QUFDZjtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsR0FBRzs7QUFFSDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7O0FBRUE7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7O0FBRUE7O0FBRUE7QUFDQTtBQUNBOztBQUVBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQSxnREFBZ0Q7QUFDaEQsK0NBQStDO0FBQy9DO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQSxXQUFXLE9BQU87QUFDbEIsV0FBVyxRQUFRO0FBQ25CLGFBQWE7QUFDYjtBQUNBOztBQUVBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTs7QUFFQTtBQUNBOztBQUVBOztBQUVBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQSxXQUFXLE9BQU87QUFDbEIsV0FBVyxPQUFPO0FBQ2xCO0FBQ0E7O0FBRUE7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQSxNQUFNO0FBQ047QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBLFdBQVcsUUFBUTtBQUNuQixXQUFXLFFBQVE7QUFDbkIsV0FBVyxRQUFRO0FBQ25CLFdBQVcsU0FBUztBQUNwQixhQUFhO0FBQ2I7QUFDQTs7QUFFQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBLFdBQVcsUUFBUTtBQUNuQixXQUFXLFFBQVE7QUFDbkIsV0FBVyxRQUFRO0FBQ25CLFdBQVcsU0FBUztBQUNwQixhQUFhO0FBQ2I7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBLDRDQUE0QztBQUM1QztBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQSx5QkFBeUI7QUFDekI7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBLFdBQVcsT0FBTztBQUNsQixhQUFhO0FBQ2I7QUFDQTs7QUFFQTs7QUFFQTtBQUNBO0FBQ0E7QUFDQTs7QUFFQTtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBOztBQUVBLElBQUk7QUFDSjtBQUNBO0FBQ0E7O0FBRUE7QUFDQTtBQUNBOztBQUVBO0FBQ0E7QUFDQTtBQUNBOztBQUVBLDhCQUE4Qjs7QUFFOUI7O0FBRUEsZ0JBQWdCO0FBQ2hCO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLElBQUk7QUFDSjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxJQUFJO0FBQ0o7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsSUFBSTtBQUNKO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFHRTs7Ozs7Ozs7Ozs7Ozs7O0FDNTRDSztBQUNQO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7Ozs7Ozs7Ozs7Ozs7Ozs7Ozs7O0FDdEk4SjtBQUN0SDtBQUNKO0FBQ3BDO0FBQ0EsdUJBQXVCLHFEQUFhLDJCQUEyQixxREFBYSwyQkFBMkIscURBQWE7QUFDcEg7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSwwRUFBMEU7QUFDMUU7QUFDQTtBQUNBLHNEQUFzRDtBQUN0RDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSwyQ0FBMkMsMENBQTBDO0FBQ3JGO0FBQ0EsdUJBQXVCO0FBQ3ZCO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSwyQkFBMkIsZ0JBQWdCO0FBQzNDO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsdUNBQXVDLDZDQUFPO0FBQzlDLHdCQUF3QixpREFBUztBQUNqQztBQUNPO0FBQ1A7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSw4Q0FBOEMsWUFBWSxrS0FBa0s7QUFDNU47QUFDQTtBQUNBO0FBQ0E7QUFDQSw4Q0FBOEMsWUFBWSxrS0FBa0s7QUFDNU47QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHFCQUFxQix1REFBZTtBQUNwQywyQkFBMkIsaUVBQXlCO0FBQ3BEO0FBQ0EsZ0RBQWdELDZDQUFPO0FBQ3ZELHdCQUF3QjtBQUN4QjtBQUNBO0FBQ0E7QUFDQSw0QkFBNEIsWUFBWTtBQUN4QztBQUNBLDBHQUEwRyxrQkFBa0I7QUFDNUgsMEdBQTBHLGtCQUFrQjtBQUM1SDtBQUNBLGlCQUFpQjtBQUNqQjtBQUNBLDhCQUE4QixrREFBVSx5RUFBeUUsOENBQVU7QUFDM0gsOEJBQThCLGtEQUFVLCtCQUErQixvREFBWSw0Q0FBNEMsa0RBQWMsY0FBYywwREFBa0I7QUFDN0s7QUFDQTtBQUNBO0FBQ0E7QUFDQSxhQUFhO0FBQ2I7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSw4QkFBOEIscURBQWE7QUFDM0M7QUFDTztBQUNQO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsOEJBQThCLHFEQUFhO0FBQzNDLG9CQUFvQiw2QkFBNkI7QUFDakQ7QUFDQTtBQUNBO0FBQ0E7QUFDTztBQUNQO0FBQ0E7QUFDQTtBQUNBOzs7Ozs7Ozs7Ozs7Ozs7QUN2S087QUFDUDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSx5RUFBeUU7QUFDekU7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSx5RUFBeUU7QUFDekU7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSx5RUFBeUU7QUFDekU7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUZBQWlGO0FBQ2pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlGQUFpRjtBQUNqRjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlGQUFpRjtBQUNqRjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlGQUFpRjtBQUNqRjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlGQUFpRjtBQUNqRjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlGQUFpRjtBQUNqRjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlGQUFpRjtBQUNqRjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlGQUFpRjtBQUNqRjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlGQUFpRjtBQUNqRjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHlGQUF5RjtBQUN6RjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlHQUFpRztBQUNqRztBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlHQUFpRztBQUNqRztBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlHQUFpRztBQUNqRztBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlGQUFpRjtBQUNqRjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlGQUFpRjtBQUNqRjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlGQUFpRjtBQUNqRjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlGQUFpRjtBQUNqRjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHlGQUF5RjtBQUN6RjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHlGQUF5RjtBQUN6RjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHlGQUF5RjtBQUN6RjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHlGQUF5RjtBQUN6RjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpRkFBaUY7QUFDakY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSx5REFBeUQ7QUFDekQ7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHlGQUF5RjtBQUN6RjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EseUZBQXlGO0FBQ3pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSx5RkFBeUY7QUFDekY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHlGQUF5RjtBQUN6RjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHlEQUF5RDtBQUN6RDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EseUZBQXlGO0FBQ3pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSx5RkFBeUY7QUFDekY7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHlGQUF5RjtBQUN6RjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EseUZBQXlGO0FBQ3pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EseURBQXlEO0FBQ3pEO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EseUZBQXlGO0FBQ3pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUVBQWlFO0FBQ2pFO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EseUZBQXlGO0FBQ3pGO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLFNBQVM7QUFDVDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxTQUFTO0FBQ1Q7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsU0FBUztBQUNUO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxTQUFTO0FBQ1Q7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxTQUFTO0FBQ1Q7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxpQkFBaUI7QUFDakI7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOzs7Ozs7Ozs7Ozs7Ozs7O0FDbDFCTztBQUNQO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNPO0FBQ1A7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOzs7Ozs7Ozs7Ozs7Ozs7Ozs7QUNqNkRvRDtBQUM3QztBQUNQO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLHVCQUF1QjtBQUN2QjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSx3QkFBd0IsVUFBVTtBQUNsQztBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLG9CQUFvQixnQkFBZ0I7QUFDcEM7QUFDQTtBQUNBO0FBQ0E7QUFDQSx3QkFBd0Isa0JBQWtCO0FBQzFDO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUNBQWlDLDRIQUE0SDtBQUM3SjtBQUNBO0FBQ0E7QUFDQSx3QkFBd0Isa0JBQWtCO0FBQzFDO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsaUNBQWlDLDRIQUE0SDtBQUM3SjtBQUNBLHNCQUFzQiwwREFBMEQ7QUFDaEY7QUFDQTtBQUNBO0FBQ0Esb0JBQW9CLHFCQUFxQjtBQUN6QztBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxTQUFTO0FBQ1Q7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLDZCQUE2QixnREFBUTtBQUNyQyxLQUFLO0FBQ0wsYUFBYTtBQUNiO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMLHVDQUF1QywrRUFBK0U7QUFDdEgsdUNBQXVDLCtFQUErRTtBQUN0SDtBQUNBO0FBQ0E7QUFDQSwwQkFBMEI7QUFDMUI7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxvQkFBb0IsdUJBQXVCO0FBQzNDO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlCQUFpQjtBQUNqQjtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsYUFBYTtBQUNiO0FBQ087QUFDUDtBQUNBO0FBQ0EsdUJBQXVCO0FBQ3ZCLHVCQUF1QjtBQUN2QjtBQUNBO0FBQ0EsUUFBUSx1REFBZTtBQUN2QixLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0EsK0JBQStCLG9CQUFvQjtBQUNuRDtBQUNPO0FBQ1A7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBLGlDQUFpQyxrQkFBa0I7QUFDbkQsNEJBQTRCO0FBQzVCO0FBQ0EsK0NBQStDO0FBQy9DO0FBQ0E7QUFDQSxLQUFLLFNBQVM7QUFDZCw0QkFBNEI7QUFDNUI7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLLFNBQVM7QUFDZDtBQUNBO0FBQ0E7QUFDQSwyQ0FBMkMsZUFBZTtBQUMxRDtBQUNBO0FBQ0Esb0JBQW9CLDRCQUE0QjtBQUNoRDtBQUNBO0FBQ0E7QUFDQTtBQUNBOzs7Ozs7Ozs7Ozs7Ozs7Ozs7Ozs7Ozs7Ozs7QUM1T3dDO0FBQ2pDO0FBQ0E7QUFDUDtBQUNBLGVBQWUsa0RBQVUsR0FBRywyREFBMkQ7QUFDdkY7QUFDTztBQUNQO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMO0FBQ087QUFDUCx1QkFBdUIsa0RBQVU7QUFDakM7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0wsYUFBYTtBQUNiO0FBQ087QUFDUDtBQUNBLG9CQUFvQixXQUFXO0FBQy9CO0FBQ0E7QUFDQTtBQUNBO0FBQ087QUFDUCwyQ0FBMkMsaUNBQWlDO0FBQzVFO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0Esb0JBQW9CLGdCQUFnQjtBQUNwQztBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDTztBQUNQO0FBQ0E7QUFDQTtBQUNBLG9CQUFvQixZQUFZO0FBQ2hDO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDTztBQUNQO0FBQ0Esb0JBQW9CLHFCQUFxQjtBQUN6Qyx3QkFBd0Isc0JBQXNCO0FBQzlDO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSx3QkFBd0Isc0JBQXNCO0FBQzlDO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSxLQUFLO0FBQ0w7QUFDTztBQUNQO0FBQ0E7QUFDQTtBQUNBO0FBQ087QUFDUDtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDTztBQUNQO0FBQ0E7QUFDTztBQUNQO0FBQ0E7QUFDQTtBQUNBO0FBQ0EsS0FBSztBQUNMOzs7Ozs7O1VDMUlBO1VBQ0E7O1VBRUE7VUFDQTtVQUNBO1VBQ0E7VUFDQTtVQUNBO1VBQ0E7VUFDQTtVQUNBO1VBQ0E7VUFDQTtVQUNBO1VBQ0E7O1VBRUE7VUFDQTs7VUFFQTtVQUNBO1VBQ0E7Ozs7O1dDdEJBO1dBQ0E7V0FDQTtXQUNBO1dBQ0EseUNBQXlDLHdDQUF3QztXQUNqRjtXQUNBO1dBQ0E7Ozs7O1dDUEE7Ozs7O1dDQUE7V0FDQTtXQUNBO1dBQ0EsdURBQXVELGlCQUFpQjtXQUN4RTtXQUNBLGdEQUFnRCxhQUFhO1dBQzdEOzs7Ozs7Ozs7Ozs7Ozs7Ozs7Ozs7OztBQ044QztBQUNhO0FBQ0o7QUFDb0MiLCJzb3VyY2VzIjpbIndlYnBhY2s6Ly9BbmltVG9Cdmgvd2VicGFjay91bml2ZXJzYWxNb2R1bGVEZWZpbml0aW9uIiwid2VicGFjazovL0FuaW1Ub0J2aC8uL25vZGVfbW9kdWxlcy9xdWF0ZXJuaW9uL2Rpc3QvcXVhdGVybmlvbi5tanMiLCJ3ZWJwYWNrOi8vQW5pbVRvQnZoLy4vc3JjL2FsaWFzZXMudHMiLCJ3ZWJwYWNrOi8vQW5pbVRvQnZoLy4vc3JjL2NvbnZlcnQudHMiLCJ3ZWJwYWNrOi8vQW5pbVRvQnZoLy4vc3JjL2hpZXJhcmNoeS50cyIsIndlYnBhY2s6Ly9BbmltVG9CdmgvLi9zcmMvb2Zmc2V0cy50cyIsIndlYnBhY2s6Ly9BbmltVG9CdmgvLi9zcmMvcGFyc2UudHMiLCJ3ZWJwYWNrOi8vQW5pbVRvQnZoLy4vc3JjL3V0aWxzLnRzIiwid2VicGFjazovL0FuaW1Ub0J2aC93ZWJwYWNrL2Jvb3RzdHJhcCIsIndlYnBhY2s6Ly9BbmltVG9Cdmgvd2VicGFjay9ydW50aW1lL2RlZmluZSBwcm9wZXJ0eSBnZXR0ZXJzIiwid2VicGFjazovL0FuaW1Ub0J2aC93ZWJwYWNrL3J1bnRpbWUvaGFzT3duUHJvcGVydHkgc2hvcnRoYW5kIiwid2VicGFjazovL0FuaW1Ub0J2aC93ZWJwYWNrL3J1bnRpbWUvbWFrZSBuYW1lc3BhY2Ugb2JqZWN0Iiwid2VicGFjazovL0FuaW1Ub0J2aC8uL3NyYy9pbmRleC50cyJdLCJzb3VyY2VzQ29udGVudCI6WyIoZnVuY3Rpb24gd2VicGFja1VuaXZlcnNhbE1vZHVsZURlZmluaXRpb24ocm9vdCwgZmFjdG9yeSkge1xuXHRpZih0eXBlb2YgZXhwb3J0cyA9PT0gJ29iamVjdCcgJiYgdHlwZW9mIG1vZHVsZSA9PT0gJ29iamVjdCcpXG5cdFx0bW9kdWxlLmV4cG9ydHMgPSBmYWN0b3J5KCk7XG5cdGVsc2UgaWYodHlwZW9mIGRlZmluZSA9PT0gJ2Z1bmN0aW9uJyAmJiBkZWZpbmUuYW1kKVxuXHRcdGRlZmluZShbXSwgZmFjdG9yeSk7XG5cdGVsc2UgaWYodHlwZW9mIGV4cG9ydHMgPT09ICdvYmplY3QnKVxuXHRcdGV4cG9ydHNbXCJBbmltVG9CdmhcIl0gPSBmYWN0b3J5KCk7XG5cdGVsc2Vcblx0XHRyb290W1wiQW5pbVRvQnZoXCJdID0gZmFjdG9yeSgpO1xufSkoc2VsZiwgKCkgPT4ge1xucmV0dXJuICIsIid1c2Ugc3RyaWN0JztcblxuLyoqXG4gKiBDcmVhdGVzIGEgbmV3IFF1YXRlcm5pb24gb2JqZWN0XG4gKiBcbiAqIEBwYXJhbSB7bnVtYmVyfSB3IFxuICogQHBhcmFtIHtudW1iZXJ9IHggXG4gKiBAcGFyYW0ge251bWJlcn0geSBcbiAqIEBwYXJhbSB7bnVtYmVyfSB6IFxuICogQHJldHVybnMgXG4gKi9cbmZ1bmN0aW9uIG5ld1F1YXRlcm5pb24odywgeCwgeSwgeikge1xuICBjb25zdCBmID0gT2JqZWN0LmNyZWF0ZShRdWF0ZXJuaW9uLnByb3RvdHlwZSk7XG5cbiAgZlsndyddID0gdztcbiAgZlsneCddID0geDtcbiAgZlsneSddID0geTtcbiAgZlsneiddID0gejtcblxuICByZXR1cm4gZjtcbn1cblxuLyoqXG4gKiBDcmVhdGVzIGEgbmV3IG5vcm1hbGl6ZWQgUXVhdGVybmlvbiBvYmplY3RcbiAqXG4gKiBAcGFyYW0ge251bWJlcn0gd1xuICogQHBhcmFtIHtudW1iZXJ9IHhcbiAqIEBwYXJhbSB7bnVtYmVyfSB5XG4gKiBAcGFyYW0ge251bWJlcn0gelxuICogQHJldHVybnNcbiAqL1xuZnVuY3Rpb24gbmV3Tm9ybWFsaXplZCh3LCB4LCB5LCB6KSB7XG4gIGNvbnN0IGYgPSBPYmplY3QuY3JlYXRlKFF1YXRlcm5pb24ucHJvdG90eXBlKTtcblxuICAvLyBXZSBhc3N1bWUgfFF8ID4gMCBmb3IgaW50ZXJuYWwgdXNhZ2VcbiAgY29uc3QgaWwgPSAxIC8gTWF0aC5zcXJ0KHcgKiB3ICsgeCAqIHggKyB5ICogeSArIHogKiB6KTtcblxuICBmWyd3J10gPSB3ICogaWw7XG4gIGZbJ3gnXSA9IHggKiBpbDtcbiAgZlsneSddID0geSAqIGlsO1xuICBmWyd6J10gPSB6ICogaWw7XG5cbiAgcmV0dXJuIGY7XG59XG5cbi8qKlxuICogQ2FsY3VsYXRlcyBsb2coc3FydChhXjIrYl4yKSkgaW4gYSB3YXkgdG8gYXZvaWQgb3ZlcmZsb3dzXG4gKlxuICogQHBhcmFtIHtudW1iZXJ9IGFcbiAqIEBwYXJhbSB7bnVtYmVyfSBiXG4gKiBAcmV0dXJucyB7bnVtYmVyfVxuICovXG5mdW5jdGlvbiBsb2dIeXBvdChhLCBiKSB7XG5cbiAgY29uc3QgX2EgPSBNYXRoLmFicyhhKTtcbiAgY29uc3QgX2IgPSBNYXRoLmFicyhiKTtcblxuICBpZiAoYSA9PT0gMCkge1xuICAgIHJldHVybiBNYXRoLmxvZyhfYik7XG4gIH1cblxuICBpZiAoYiA9PT0gMCkge1xuICAgIHJldHVybiBNYXRoLmxvZyhfYSk7XG4gIH1cblxuICBpZiAoX2EgPCAzMDAwICYmIF9iIDwgMzAwMCkge1xuICAgIHJldHVybiAwLjUgKiBNYXRoLmxvZyhhICogYSArIGIgKiBiKTtcbiAgfVxuXG4gIGEgPSBhIC8gMjtcbiAgYiA9IGIgLyAyO1xuXG4gIHJldHVybiAwLjUgKiBNYXRoLmxvZyhhICogYSArIGIgKiBiKSArIE1hdGguTE4yO1xufVxuXG4vKlxuICogVGVtcG9yYXJ5IHBhcnNpbmcgb2JqZWN0IHRvIGF2b2lkIHJlLWFsbG9jYXRpb25zXG4gKlxuICovXG5jb25zdCBQID0gT2JqZWN0LmNyZWF0ZShRdWF0ZXJuaW9uLnByb3RvdHlwZSk7XG5cbmZ1bmN0aW9uIHBhcnNlKGRlc3QsIHcsIHgsIHksIHopIHtcblxuICAvLyBNb3N0IGNvbW1vbiBpbnRlcm5hbCB1c2UgY2FzZSB3aXRoIDQgcGFyYW1zXG4gIGlmICh6ICE9PSB1bmRlZmluZWQpIHtcbiAgICBkZXN0Wyd3J10gPSB3O1xuICAgIGRlc3RbJ3gnXSA9IHg7XG4gICAgZGVzdFsneSddID0geTtcbiAgICBkZXN0Wyd6J10gPSB6O1xuICAgIHJldHVybjtcbiAgfVxuXG4gIGlmICh0eXBlb2YgdyA9PT0gJ29iamVjdCcgJiYgeSA9PT0gdW5kZWZpbmVkKSB7XG5cbiAgICAvLyBDaGVjayBmb3IgcXVhdHMsIGZvciBleGFtcGxlIHdoZW4gYW4gb2JqZWN0IGdldHMgY2xvbmVkXG4gICAgaWYgKCd3JyBpbiB3IHx8ICd4JyBpbiB3IHx8ICd5JyBpbiB3IHx8ICd6JyBpbiB3KSB7XG4gICAgICBkZXN0Wyd3J10gPSB3Wyd3J10gfHwgMDtcbiAgICAgIGRlc3RbJ3gnXSA9IHdbJ3gnXSB8fCAwO1xuICAgICAgZGVzdFsneSddID0gd1sneSddIHx8IDA7XG4gICAgICBkZXN0Wyd6J10gPSB3Wyd6J10gfHwgMDtcbiAgICAgIHJldHVybjtcbiAgICB9XG5cbiAgICAvLyBDaGVjayBmb3IgY29tcGxleCBudW1iZXJzXG4gICAgaWYgKCdyZScgaW4gdyAmJiAnaW0nIGluIHcpIHtcbiAgICAgIGRlc3RbJ3cnXSA9IHdbJ3JlJ107XG4gICAgICBkZXN0Wyd4J10gPSB3WydpbSddO1xuICAgICAgZGVzdFsneSddID0gMDtcbiAgICAgIGRlc3RbJ3onXSA9IDA7XG4gICAgICByZXR1cm47XG4gICAgfVxuXG4gICAgLy8gQ2hlY2sgZm9yIGFycmF5XG4gICAgaWYgKHcubGVuZ3RoID09PSA0KSB7XG4gICAgICBkZXN0Wyd3J10gPSB3WzBdO1xuICAgICAgZGVzdFsneCddID0gd1sxXTtcbiAgICAgIGRlc3RbJ3knXSA9IHdbMl07XG4gICAgICBkZXN0Wyd6J10gPSB3WzNdO1xuICAgICAgcmV0dXJuO1xuICAgIH1cblxuICAgIC8vIENoZWNrIGZvciBhdWdtZW50ZWQgdmVjdG9yXG4gICAgaWYgKHcubGVuZ3RoID09PSAzKSB7XG4gICAgICBkZXN0Wyd3J10gPSAwO1xuICAgICAgZGVzdFsneCddID0gd1swXTtcbiAgICAgIGRlc3RbJ3knXSA9IHdbMV07XG4gICAgICBkZXN0Wyd6J10gPSB3WzJdO1xuICAgICAgcmV0dXJuO1xuICAgIH1cblxuICAgIHRocm93IG5ldyBFcnJvcignSW52YWxpZCBvYmplY3QnKTtcbiAgfVxuXG4gIC8vIFBhcnNlIHN0cmluZyB2YWx1ZXNcbiAgaWYgKHR5cGVvZiB3ID09PSAnc3RyaW5nJyAmJiB5ID09PSB1bmRlZmluZWQpIHtcblxuICAgIGNvbnN0IHRva2VucyA9IHcubWF0Y2goL1xcZCtcXC4/XFxkKmVbKy1dP1xcZCt8XFxkK1xcLj9cXGQqfFxcLlxcZCt8Li9nKTtcbiAgICBsZXQgcGx1cyA9IDE7XG4gICAgbGV0IG1pbnVzID0gMDtcblxuICAgIGNvbnN0IGlNYXAgPSB7ICdpJzogJ3gnLCAnaic6ICd5JywgJ2snOiAneicgfTtcblxuICAgIGlmICh0b2tlbnMgPT09IG51bGwpIHtcbiAgICAgIHRocm93IG5ldyBFcnJvcignUGFyc2UgZXJyb3InKTtcbiAgICB9XG5cbiAgICAvLyBSZXNldCB0aGUgY3VycmVudCBzdGF0ZVxuICAgIGRlc3RbJ3cnXSA9XG4gICAgICBkZXN0Wyd4J10gPVxuICAgICAgZGVzdFsneSddID1cbiAgICAgIGRlc3RbJ3onXSA9IDA7XG5cbiAgICBmb3IgKGxldCBpID0gMDsgaSA8IHRva2Vucy5sZW5ndGg7IGkrKykge1xuXG4gICAgICBsZXQgYyA9IHRva2Vuc1tpXTtcbiAgICAgIGxldCBkID0gdG9rZW5zW2kgKyAxXTtcblxuICAgICAgaWYgKGMgPT09ICcgJyB8fCBjID09PSAnXFx0JyB8fCBjID09PSAnXFxuJykge1xuICAgICAgICAvKiB2b2lkICovXG4gICAgICB9IGVsc2UgaWYgKGMgPT09ICcrJykge1xuICAgICAgICBwbHVzKys7XG4gICAgICB9IGVsc2UgaWYgKGMgPT09ICctJykge1xuICAgICAgICBtaW51cysrO1xuICAgICAgfSBlbHNlIHtcblxuICAgICAgICBpZiAocGx1cyArIG1pbnVzID09PSAwKSB7XG4gICAgICAgICAgdGhyb3cgbmV3IEVycm9yKCdQYXJzZSBlcnJvcicgKyBjKTtcbiAgICAgICAgfVxuICAgICAgICBsZXQgZyA9IGlNYXBbY107XG5cbiAgICAgICAgLy8gSXMgdGhlIGN1cnJlbnQgdG9rZW4gYW4gaW1hZ2luYXJ5IHNpZ24/XG4gICAgICAgIGlmIChnICE9PSB1bmRlZmluZWQpIHtcblxuICAgICAgICAgIC8vIElzIHRoZSBmb2xsb3dpbmcgdG9rZW4gYSBudW1iZXI/XG4gICAgICAgICAgaWYgKGQgIT09ICcgJyAmJiAhaXNOYU4oZCkpIHtcbiAgICAgICAgICAgIGMgPSBkO1xuICAgICAgICAgICAgaSsrO1xuICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBjID0gJzEnO1xuICAgICAgICAgIH1cblxuICAgICAgICB9IGVsc2Uge1xuXG4gICAgICAgICAgaWYgKGlzTmFOKGMpKSB7XG4gICAgICAgICAgICB0aHJvdyBuZXcgRXJyb3IoJ1BhcnNlciBlcnJvcicpO1xuICAgICAgICAgIH1cblxuICAgICAgICAgIGcgPSBpTWFwW2RdO1xuXG4gICAgICAgICAgaWYgKGcgIT09IHVuZGVmaW5lZCkge1xuICAgICAgICAgICAgaSsrO1xuICAgICAgICAgIH1cbiAgICAgICAgfVxuXG4gICAgICAgIGRlc3RbZyB8fCAndyddICs9IHBhcnNlRmxvYXQoKG1pbnVzICUgMiA/ICctJyA6ICcnKSArIGMpO1xuICAgICAgICBwbHVzID0gbWludXMgPSAwO1xuICAgICAgfVxuICAgIH1cblxuICAgIC8vIFN0aWxsIHNvbWV0aGluZyBvbiB0aGUgc3RhY2tcbiAgICBpZiAocGx1cyArIG1pbnVzID4gMCkge1xuICAgICAgdGhyb3cgbmV3IEVycm9yKCdQYXJzZXIgZXJyb3InKTtcbiAgICB9XG4gICAgcmV0dXJuO1xuICB9XG5cbiAgLy8gSWYgbm8gc2luZ2xlIHZhcmlhYmxlIHdhcyBnaXZlbiBBTkQgaXQgd2FzIHRoZSBjb25zdHJ1Y3Rvciwgc2V0IGl0IHRvIHRoZSBpZGVudGl0eVxuICBpZiAodyA9PT0gdW5kZWZpbmVkICYmIGRlc3QgIT09IFApIHtcbiAgICBkZXN0Wyd3J10gPSAxO1xuICAgIGRlc3RbJ3gnXSA9XG4gICAgICBkZXN0Wyd5J10gPVxuICAgICAgZGVzdFsneiddID0gMDtcbiAgfSBlbHNlIHtcblxuICAgIGRlc3RbJ3cnXSA9IHcgfHwgMDtcblxuICAgIC8vIE5vdGU6IFRoaXMgaXNuJ3QgZnJvbUF4aXMoKSwgaXQncyBqdXN0IHN5bnRhY3RpYyBzdWdhciFcbiAgICBpZiAoeCAmJiB4Lmxlbmd0aCA9PT0gMykge1xuICAgICAgZGVzdFsneCddID0geFswXTtcbiAgICAgIGRlc3RbJ3knXSA9IHhbMV07XG4gICAgICBkZXN0Wyd6J10gPSB4WzJdO1xuICAgIH0gZWxzZSB7XG4gICAgICBkZXN0Wyd4J10gPSB4IHx8IDA7XG4gICAgICBkZXN0Wyd5J10gPSB5IHx8IDA7XG4gICAgICBkZXN0Wyd6J10gPSB6IHx8IDA7XG4gICAgfVxuICB9XG59XG5cbmZ1bmN0aW9uIG51bVRvU3RyKG4sIGNoYXIsIHByZXYpIHtcblxuICBsZXQgcmV0ID0gJyc7XG5cbiAgaWYgKG4gIT09IDApIHtcblxuICAgIGlmIChwcmV2ICE9PSAnJykge1xuICAgICAgcmV0ICs9IG4gPCAwID8gJyAtICcgOiAnICsgJztcbiAgICB9IGVsc2UgaWYgKG4gPCAwKSB7XG4gICAgICByZXQgKz0gJy0nO1xuICAgIH1cblxuICAgIG4gPSBNYXRoLmFicyhuKTtcblxuICAgIGlmICgxICE9PSBuIHx8IGNoYXIgPT09ICcnKSB7XG4gICAgICByZXQgKz0gbjtcbiAgICB9XG4gICAgcmV0ICs9IGNoYXI7XG4gIH1cbiAgcmV0dXJuIHJldDtcbn1cblxuLyoqXG4gKiBRdWF0ZXJuaW9uIGNvbnN0cnVjdG9yXG4gKlxuICogQGNvbnN0cnVjdG9yXG4gKiBAcGFyYW0ge251bWJlcnxPYmplY3R8c3RyaW5nfSB3IHJlYWxcbiAqIEBwYXJhbSB7bnVtYmVyPX0geCBpbWFnXG4gKiBAcGFyYW0ge251bWJlcj19IHkgaW1hZ1xuICogQHBhcmFtIHtudW1iZXI9fSB6IGltYWdcbiAqIEByZXR1cm5zIHtRdWF0ZXJuaW9ufVxuICovXG5mdW5jdGlvbiBRdWF0ZXJuaW9uKHcsIHgsIHksIHopIHtcblxuICBpZiAodGhpcyBpbnN0YW5jZW9mIFF1YXRlcm5pb24pIHtcbiAgICBwYXJzZSh0aGlzLCB3LCB4LCB5LCB6KTtcbiAgfSBlbHNlIHtcbiAgICBjb25zdCB0ID0gT2JqZWN0LmNyZWF0ZShRdWF0ZXJuaW9uLnByb3RvdHlwZSk7XG4gICAgcGFyc2UodCwgdywgeCwgeSwgeik7XG4gICAgcmV0dXJuIHQ7XG4gIH1cbn1cblxuUXVhdGVybmlvbi5wcm90b3R5cGUgPSB7XG4gICd3JzogMSxcbiAgJ3gnOiAwLFxuICAneSc6IDAsXG4gICd6JzogMCxcbiAgLyoqXG4gICAqIEFkZHMgdHdvIHF1YXRlcm5pb25zIFExIGFuZCBRMlxuICAgKlxuICAgKiBAcGFyYW0ge251bWJlcnxPYmplY3R8c3RyaW5nfSB3IHJlYWxcbiAgICogQHBhcmFtIHtudW1iZXI9fSB4IGltYWdcbiAgICogQHBhcmFtIHtudW1iZXI9fSB5IGltYWdcbiAgICogQHBhcmFtIHtudW1iZXI9fSB6IGltYWdcbiAgICogQHJldHVybnMge1F1YXRlcm5pb259XG4gICAqL1xuICAnYWRkJzogZnVuY3Rpb24gKHcsIHgsIHksIHopIHtcblxuICAgIHBhcnNlKFAsIHcsIHgsIHksIHopO1xuXG4gICAgLy8gUTEgKyBRMiA6PSBbdzEsIHYxXSArIFt3MiwgdjJdID0gW3cxICsgdzIsIHYxICsgdjJdXG5cbiAgICByZXR1cm4gbmV3UXVhdGVybmlvbihcbiAgICAgIHRoaXNbJ3cnXSArIFBbJ3cnXSxcbiAgICAgIHRoaXNbJ3gnXSArIFBbJ3gnXSxcbiAgICAgIHRoaXNbJ3knXSArIFBbJ3knXSxcbiAgICAgIHRoaXNbJ3onXSArIFBbJ3onXSk7XG4gIH0sXG4gIC8qKlxuICAgKiBTdWJ0cmFjdHMgYSBxdWF0ZXJuaW9ucyBRMiBmcm9tIFExXG4gICAqXG4gICAqIEBwYXJhbSB7bnVtYmVyfE9iamVjdHxzdHJpbmd9IHcgcmVhbFxuICAgKiBAcGFyYW0ge251bWJlcj19IHggaW1hZ1xuICAgKiBAcGFyYW0ge251bWJlcj19IHkgaW1hZ1xuICAgKiBAcGFyYW0ge251bWJlcj19IHogaW1hZ1xuICAgKiBAcmV0dXJucyB7UXVhdGVybmlvbn1cbiAgICovXG4gICdzdWInOiBmdW5jdGlvbiAodywgeCwgeSwgeikge1xuXG4gICAgcGFyc2UoUCwgdywgeCwgeSwgeik7XG5cbiAgICAvLyBRMSAtIFEyIDo9IFExICsgKC1RMilcbiAgICAvLyAgICAgICAgICA9IFt3MSwgdjFdIC0gW3cyLCB2Ml0gPSBbdzEgLSB3MiwgdjEgLSB2Ml1cblxuICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKFxuICAgICAgdGhpc1sndyddIC0gUFsndyddLFxuICAgICAgdGhpc1sneCddIC0gUFsneCddLFxuICAgICAgdGhpc1sneSddIC0gUFsneSddLFxuICAgICAgdGhpc1sneiddIC0gUFsneiddKTtcbiAgfSxcbiAgLyoqXG4gICAqIENhbGN1bGF0ZXMgdGhlIGFkZGl0aXZlIGludmVyc2UsIG9yIHNpbXBseSBpdCBuZWdhdGVzIHRoZSBxdWF0ZXJuaW9uXG4gICAqXG4gICAqIEByZXR1cm5zIHtRdWF0ZXJuaW9ufVxuICAgKi9cbiAgJ25lZyc6IGZ1bmN0aW9uICgpIHtcblxuICAgIC8vIC1RIDo9IFstdywgLXZdXG5cbiAgICByZXR1cm4gbmV3UXVhdGVybmlvbigtdGhpc1sndyddLCAtdGhpc1sneCddLCAtdGhpc1sneSddLCAtdGhpc1sneiddKTtcbiAgfSxcbiAgLyoqXG4gICAqIENhbGN1bGF0ZXMgdGhlIGxlbmd0aC9tb2R1bHVzL21hZ25pdHVkZSBvciB0aGUgbm9ybSBvZiBhIHF1YXRlcm5pb25cbiAgICpcbiAgICogQHJldHVybnMge251bWJlcn1cbiAgICovXG4gICdub3JtJzogZnVuY3Rpb24gKCkge1xuXG4gICAgLy8gfFF8IDo9IHNxcnQofFF8XjIpXG5cbiAgICAvLyBUaGUgdW5pdCBxdWF0ZXJuaW9uIGhhcyB8UXwgPSAxXG5cbiAgICBjb25zdCB3ID0gdGhpc1sndyddO1xuICAgIGNvbnN0IHggPSB0aGlzWyd4J107XG4gICAgY29uc3QgeSA9IHRoaXNbJ3knXTtcbiAgICBjb25zdCB6ID0gdGhpc1sneiddO1xuXG4gICAgcmV0dXJuIE1hdGguc3FydCh3ICogdyArIHggKiB4ICsgeSAqIHkgKyB6ICogeik7XG4gIH0sXG4gIC8qKlxuICAgKiBDYWxjdWxhdGVzIHRoZSBzcXVhcmVkIGxlbmd0aC9tb2R1bHVzL21hZ25pdHVkZSBvciB0aGUgbm9ybSBvZiBhIHF1YXRlcm5pb25cbiAgICpcbiAgICogQHJldHVybnMge251bWJlcn1cbiAgICovXG4gICdub3JtU3EnOiBmdW5jdGlvbiAoKSB7XG5cbiAgICAvLyB8UXxeMiA6PSBbdywgdl0gKiBbdywgLXZdXG4gICAgLy8gICAgICAgID0gW3deMiArIGRvdCh2LCB2KSwgLXcgKiB2ICsgdyAqIHYgKyBjcm9zcyh2LCAtdildXG4gICAgLy8gICAgICAgID0gW3deMiArIHx2fF4yLCAwXVxuICAgIC8vICAgICAgICA9IFt3XjIgKyBkb3QodiwgdiksIDBdXG4gICAgLy8gICAgICAgID0gZG90KFEsIFEpXG4gICAgLy8gICAgICAgID0gUSAqIFEnXG5cbiAgICBjb25zdCB3ID0gdGhpc1sndyddO1xuICAgIGNvbnN0IHggPSB0aGlzWyd4J107XG4gICAgY29uc3QgeSA9IHRoaXNbJ3knXTtcbiAgICBjb25zdCB6ID0gdGhpc1sneiddO1xuXG4gICAgcmV0dXJuIHcgKiB3ICsgeCAqIHggKyB5ICogeSArIHogKiB6O1xuICB9LFxuICAvKipcbiAgICogTm9ybWFsaXplcyB0aGUgcXVhdGVybmlvbiB0byBoYXZlIHxRfCA9IDEgYXMgbG9uZyBhcyB0aGUgbm9ybSBpcyBub3QgemVyb1xuICAgKiBBbHRlcm5hdGl2ZSBuYW1lcyBhcmUgdGhlIHNpZ251bSwgdW5pdCBvciB2ZXJzb3JcbiAgICpcbiAgICogQHJldHVybnMge1F1YXRlcm5pb259XG4gICAqL1xuICAnbm9ybWFsaXplJzogZnVuY3Rpb24gKCkge1xuXG4gICAgLy8gUSogOj0gUSAvIHxRfFxuXG4gICAgLy8gdW5yb2xsZWQgUS5zY2FsZSgxIC8gUS5ub3JtKCkpXG5cbiAgICBjb25zdCB3ID0gdGhpc1sndyddO1xuICAgIGNvbnN0IHggPSB0aGlzWyd4J107XG4gICAgY29uc3QgeSA9IHRoaXNbJ3knXTtcbiAgICBjb25zdCB6ID0gdGhpc1sneiddO1xuXG4gICAgbGV0IG5vcm0gPSBNYXRoLnNxcnQodyAqIHcgKyB4ICogeCArIHkgKiB5ICsgeiAqIHopO1xuXG4gICAgaWYgKG5vcm0gPCBFUFNJTE9OKSB7XG4gICAgICByZXR1cm4gUXVhdGVybmlvblsnWkVSTyddO1xuICAgIH1cblxuICAgIG5vcm0gPSAxIC8gbm9ybTtcblxuICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKHcgKiBub3JtLCB4ICogbm9ybSwgeSAqIG5vcm0sIHogKiBub3JtKTtcbiAgfSxcbiAgLyoqXG4gICAqIENhbGN1bGF0ZXMgdGhlIEhhbWlsdG9uIHByb2R1Y3Qgb2YgdHdvIHF1YXRlcm5pb25zXG4gICAqIExlYXZpbmcgb3V0IHRoZSBpbWFnaW5hcnkgcGFydCByZXN1bHRzIGluIGp1c3Qgc2NhbGluZyB0aGUgcXVhdFxuICAgKlxuICAgKiBAcGFyYW0ge251bWJlcnxPYmplY3R8c3RyaW5nfSB3IHJlYWxcbiAgICogQHBhcmFtIHtudW1iZXI9fSB4IGltYWdcbiAgICogQHBhcmFtIHtudW1iZXI9fSB5IGltYWdcbiAgICogQHBhcmFtIHtudW1iZXI9fSB6IGltYWdcbiAgICogQHJldHVybnMge1F1YXRlcm5pb259XG4gICAqL1xuICAnbXVsJzogZnVuY3Rpb24gKHcsIHgsIHksIHopIHtcblxuICAgIHBhcnNlKFAsIHcsIHgsIHksIHopO1xuXG4gICAgLy8gUTEgKiBRMiA9IFt3MSAqIHcyIC0gZG90KHYxLCB2MiksIHcxICogdjIgKyB3MiAqIHYxICsgY3Jvc3ModjEsIHYyKV1cblxuICAgIC8vIE5vdCBjb21tdXRhdGl2ZSBiZWNhdXNlIGNyb3NzKHYxLCB2MikgIT0gY3Jvc3ModjIsIHYxKSFcblxuICAgIGNvbnN0IHcxID0gdGhpc1sndyddO1xuICAgIGNvbnN0IHgxID0gdGhpc1sneCddO1xuICAgIGNvbnN0IHkxID0gdGhpc1sneSddO1xuICAgIGNvbnN0IHoxID0gdGhpc1sneiddO1xuXG4gICAgY29uc3QgdzIgPSBQWyd3J107XG4gICAgY29uc3QgeDIgPSBQWyd4J107XG4gICAgY29uc3QgeTIgPSBQWyd5J107XG4gICAgY29uc3QgejIgPSBQWyd6J107XG5cbiAgICByZXR1cm4gbmV3UXVhdGVybmlvbihcbiAgICAgIHcxICogdzIgLSB4MSAqIHgyIC0geTEgKiB5MiAtIHoxICogejIsXG4gICAgICB3MSAqIHgyICsgeDEgKiB3MiArIHkxICogejIgLSB6MSAqIHkyLFxuICAgICAgdzEgKiB5MiArIHkxICogdzIgKyB6MSAqIHgyIC0geDEgKiB6MixcbiAgICAgIHcxICogejIgKyB6MSAqIHcyICsgeDEgKiB5MiAtIHkxICogeDIpO1xuICB9LFxuICAvKipcbiAgICogU2NhbGVzIGEgcXVhdGVybmlvbiBieSBhIHNjYWxhciwgZmFzdGVyIHRoYW4gdXNpbmcgbXVsdGlwbGljYXRpb25cbiAgICpcbiAgICogQHBhcmFtIHtudW1iZXJ9IHMgc2NhbGluZyBmYWN0b3JcbiAgICogQHJldHVybnMge1F1YXRlcm5pb259XG4gICAqL1xuICAnc2NhbGUnOiBmdW5jdGlvbiAocykge1xuXG4gICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oXG4gICAgICB0aGlzWyd3J10gKiBzLFxuICAgICAgdGhpc1sneCddICogcyxcbiAgICAgIHRoaXNbJ3knXSAqIHMsXG4gICAgICB0aGlzWyd6J10gKiBzKTtcbiAgfSxcbiAgLyoqXG4gICAqIENhbGN1bGF0ZXMgdGhlIGRvdCBwcm9kdWN0IG9mIHR3byBxdWF0ZXJuaW9uc1xuICAgKlxuICAgKiBAcGFyYW0ge251bWJlcnxPYmplY3R8c3RyaW5nfSB3IHJlYWxcbiAgICogQHBhcmFtIHtudW1iZXI9fSB4IGltYWdcbiAgICogQHBhcmFtIHtudW1iZXI9fSB5IGltYWdcbiAgICogQHBhcmFtIHtudW1iZXI9fSB6IGltYWdcbiAgICogQHJldHVybnMge251bWJlcn1cbiAgICovXG4gICdkb3QnOiBmdW5jdGlvbiAodywgeCwgeSwgeikge1xuXG4gICAgcGFyc2UoUCwgdywgeCwgeSwgeik7XG5cbiAgICAvLyBkb3QoUTEsIFEyKSA6PSB3MSAqIHcyICsgZG90KHYxLCB2MilcblxuICAgIHJldHVybiB0aGlzWyd3J10gKiBQWyd3J10gKyB0aGlzWyd4J10gKiBQWyd4J10gKyB0aGlzWyd5J10gKiBQWyd5J10gKyB0aGlzWyd6J10gKiBQWyd6J107XG4gIH0sXG4gIC8qKlxuICAgKiBDYWxjdWxhdGVzIHRoZSBpbnZlcnNlIG9mIGEgcXVhdCBmb3Igbm9uLW5vcm1hbGl6ZWQgcXVhdHMgc3VjaCB0aGF0XG4gICAqIFFeLTEgKiBRID0gMSBhbmQgUSAqIFFeLTEgPSAxXG4gICAqXG4gICAqIEByZXR1cm5zIHtRdWF0ZXJuaW9ufVxuICAgKi9cbiAgJ2ludmVyc2UnOiBmdW5jdGlvbiAoKSB7XG5cbiAgICAvLyBRXi0xIDo9IFEnIC8gfFF8XjJcbiAgICAvLyAgICAgICA9IFt3IC8gKHdeMiArIHx2fF4yKSwgLXYgLyAod14yICsgfHZ8XjIpXVxuXG4gICAgLy8gUHJvb2Y6XG4gICAgLy8gUSAqIFFeLTEgPSBbdywgdl0gKiBbdyAvICh3XjIgKyB8dnxeMiksIC12IC8gKHdeMiArIHx2fF4yKV1cbiAgICAvLyAgICAgICAgICA9IFsxLCAwXVxuICAgIC8vIFFeLTEgKiBRID0gW3cgLyAod14yICsgfHZ8XjIpLCAtdiAvICh3XjIgKyB8dnxeMildICogW3csIHZdXG4gICAgLy8gICAgICAgICAgPSBbMSwgMF0uXG5cbiAgICBjb25zdCB3ID0gdGhpc1sndyddO1xuICAgIGNvbnN0IHggPSB0aGlzWyd4J107XG4gICAgY29uc3QgeSA9IHRoaXNbJ3knXTtcbiAgICBjb25zdCB6ID0gdGhpc1sneiddO1xuXG4gICAgbGV0IG5vcm1TcSA9IHcgKiB3ICsgeCAqIHggKyB5ICogeSArIHogKiB6O1xuXG4gICAgaWYgKG5vcm1TcSA9PT0gMCkge1xuICAgICAgcmV0dXJuIFF1YXRlcm5pb25bJ1pFUk8nXTsgLy8gVE9ETzogSXMgdGhlIHJlc3VsdCB6ZXJvIG9yIG9uZSB3aGVuIHRoZSBub3JtIGlzIHplcm8/XG4gICAgfVxuXG4gICAgbm9ybVNxID0gMSAvIG5vcm1TcTtcblxuICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKHcgKiBub3JtU3EsIC14ICogbm9ybVNxLCAteSAqIG5vcm1TcSwgLXogKiBub3JtU3EpO1xuICB9LFxuICAvKipcbiAgICogTXVsdGlwbGllcyBhIHF1YXRlcm5pb24gd2l0aCB0aGUgaW52ZXJzZSBvZiBhIHNlY29uZCBxdWF0ZXJuaW9uXG4gICAqXG4gICAqIEBwYXJhbSB7bnVtYmVyfE9iamVjdHxzdHJpbmd9IHcgcmVhbFxuICAgKiBAcGFyYW0ge251bWJlcj19IHggaW1hZ1xuICAgKiBAcGFyYW0ge251bWJlcj19IHkgaW1hZ1xuICAgKiBAcGFyYW0ge251bWJlcj19IHogaW1hZ1xuICAgKiBAcmV0dXJucyB7UXVhdGVybmlvbn1cbiAgICovXG4gICdkaXYnOiBmdW5jdGlvbiAodywgeCwgeSwgeikge1xuXG4gICAgcGFyc2UoUCwgdywgeCwgeSwgeik7XG5cbiAgICAvLyBRMSAvIFEyIDo9IFExICogUTJeLTFcblxuICAgIGNvbnN0IHcxID0gdGhpc1sndyddO1xuICAgIGNvbnN0IHgxID0gdGhpc1sneCddO1xuICAgIGNvbnN0IHkxID0gdGhpc1sneSddO1xuICAgIGNvbnN0IHoxID0gdGhpc1sneiddO1xuXG4gICAgY29uc3QgdzIgPSBQWyd3J107XG4gICAgY29uc3QgeDIgPSBQWyd4J107XG4gICAgY29uc3QgeTIgPSBQWyd5J107XG4gICAgY29uc3QgejIgPSBQWyd6J107XG5cbiAgICBsZXQgbm9ybVNxID0gdzIgKiB3MiArIHgyICogeDIgKyB5MiAqIHkyICsgejIgKiB6MjtcblxuICAgIGlmIChub3JtU3EgPT09IDApIHtcbiAgICAgIHJldHVybiBRdWF0ZXJuaW9uWydaRVJPJ107IC8vIFRPRE86IElzIHRoZSByZXN1bHQgemVybyBvciBvbmUgd2hlbiB0aGUgbm9ybSBpcyB6ZXJvP1xuICAgIH1cblxuICAgIG5vcm1TcSA9IDEgLyBub3JtU3E7XG5cbiAgICByZXR1cm4gbmV3UXVhdGVybmlvbihcbiAgICAgICh3MSAqIHcyICsgeDEgKiB4MiArIHkxICogeTIgKyB6MSAqIHoyKSAqIG5vcm1TcSxcbiAgICAgICh4MSAqIHcyIC0gdzEgKiB4MiAtIHkxICogejIgKyB6MSAqIHkyKSAqIG5vcm1TcSxcbiAgICAgICh5MSAqIHcyIC0gdzEgKiB5MiAtIHoxICogeDIgKyB4MSAqIHoyKSAqIG5vcm1TcSxcbiAgICAgICh6MSAqIHcyIC0gdzEgKiB6MiAtIHgxICogeTIgKyB5MSAqIHgyKSAqIG5vcm1TcSk7XG4gIH0sXG4gIC8qKlxuICAgKiBDYWxjdWxhdGVzIHRoZSBjb25qdWdhdGUgb2YgYSBxdWF0ZXJuaW9uXG4gICAqXG4gICAqIEByZXR1cm5zIHtRdWF0ZXJuaW9ufVxuICAgKi9cbiAgJ2Nvbmp1Z2F0ZSc6IGZ1bmN0aW9uICgpIHtcblxuICAgIC8vIFEnID0gW3MsIC12XVxuXG4gICAgLy8gSWYgdGhlIHF1YXRlcm5pb24gaXMgbm9ybWFsaXplZCxcbiAgICAvLyB0aGUgY29uanVnYXRlIGlzIHRoZSBpbnZlcnNlIG9mIHRoZSBxdWF0ZXJuaW9uIC0gYnV0IGZhc3RlclxuICAgIC8vIFEnICogUSA9IFEgKiBRJyA9IDFcblxuICAgIC8vIEFkZGl0aW9uYWxseSwgdGhlIGNvbmp1Z2F0ZSBvZiBhIHVuaXQgcXVhdGVybmlvbiBpcyBhIHJvdGF0aW9uIHdpdGggdGhlIHNhbWVcbiAgICAvLyBhbmdsZSBidXQgdGhlIG9wcG9zaXRlIGF4aXMuXG5cbiAgICAvLyBNb3Jlb3ZlciB0aGUgZm9sbG93aW5nIHByb3BlcnR5IGhvbGRzOlxuICAgIC8vIChRMSAqIFEyKScgPSBRMicgKiBRMSdcblxuICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKHRoaXNbJ3cnXSwgLXRoaXNbJ3gnXSwgLXRoaXNbJ3knXSwgLXRoaXNbJ3onXSk7XG4gIH0sXG4gIC8qKlxuICAgKiBDYWxjdWxhdGVzIHRoZSBuYXR1cmFsIGV4cG9uZW50aWF0aW9uIG9mIHRoZSBxdWF0ZXJuaW9uXG4gICAqXG4gICAqIEByZXR1cm5zIHtRdWF0ZXJuaW9ufVxuICAgKi9cbiAgJ2V4cCc6IGZ1bmN0aW9uICgpIHtcblxuICAgIGNvbnN0IHcgPSB0aGlzWyd3J107XG4gICAgY29uc3QgeCA9IHRoaXNbJ3gnXTtcbiAgICBjb25zdCB5ID0gdGhpc1sneSddO1xuICAgIGNvbnN0IHogPSB0aGlzWyd6J107XG5cbiAgICBjb25zdCB2Tm9ybSA9IE1hdGguc3FydCh4ICogeCArIHkgKiB5ICsgeiAqIHopO1xuICAgIGNvbnN0IHdFeHAgPSBNYXRoLmV4cCh3KTtcbiAgICBjb25zdCBzY2FsZSA9IHdFeHAgKiBNYXRoLnNpbih2Tm9ybSkgLyB2Tm9ybTtcblxuICAgIGlmICh2Tm9ybSA9PT0gMCkge1xuICAgICAgLy9yZXR1cm4gbmV3UXVhdGVybmlvbih3RXhwICogTWF0aC5jb3Modk5vcm0pLCAwLCAwLCAwKTtcbiAgICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKHdFeHAsIDAsIDAsIDApO1xuICAgIH1cblxuICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKFxuICAgICAgd0V4cCAqIE1hdGguY29zKHZOb3JtKSxcbiAgICAgIHggKiBzY2FsZSxcbiAgICAgIHkgKiBzY2FsZSxcbiAgICAgIHogKiBzY2FsZSk7XG4gIH0sXG4gIC8qKlxuICAgKiBDYWxjdWxhdGVzIHRoZSBuYXR1cmFsIGxvZ2FyaXRobSBvZiB0aGUgcXVhdGVybmlvblxuICAgKlxuICAgKiBAcmV0dXJucyB7UXVhdGVybmlvbn1cbiAgICovXG4gICdsb2cnOiBmdW5jdGlvbiAoKSB7XG5cbiAgICBjb25zdCB3ID0gdGhpc1sndyddO1xuICAgIGNvbnN0IHggPSB0aGlzWyd4J107XG4gICAgY29uc3QgeSA9IHRoaXNbJ3knXTtcbiAgICBjb25zdCB6ID0gdGhpc1sneiddO1xuXG4gICAgaWYgKHkgPT09IDAgJiYgeiA9PT0gMCkge1xuICAgICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oXG4gICAgICAgIGxvZ0h5cG90KHcsIHgpLFxuICAgICAgICBNYXRoLmF0YW4yKHgsIHcpLCAwLCAwKTtcbiAgICB9XG5cbiAgICBjb25zdCBxTm9ybTIgPSB4ICogeCArIHkgKiB5ICsgeiAqIHogKyB3ICogdztcbiAgICBjb25zdCB2Tm9ybSA9IE1hdGguc3FydCh4ICogeCArIHkgKiB5ICsgeiAqIHopO1xuXG4gICAgY29uc3Qgc2NhbGUgPSBNYXRoLmF0YW4yKHZOb3JtLCB3KSAvIHZOb3JtOyAvLyBBbHRlcm5hdGl2ZTogYWNvcyh3IC8gcU5vcm0pIC8gdk5vcm1cblxuICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKFxuICAgICAgTWF0aC5sb2cocU5vcm0yKSAqIDAuNSxcbiAgICAgIHggKiBzY2FsZSxcbiAgICAgIHkgKiBzY2FsZSxcbiAgICAgIHogKiBzY2FsZSk7XG4gIH0sXG4gIC8qKlxuICAgKiBDYWxjdWxhdGVzIHRoZSBwb3dlciBvZiBhIHF1YXRlcm5pb24gcmFpc2VkIHRvIGEgcmVhbCBudW1iZXIgb3IgYW5vdGhlciBxdWF0ZXJuaW9uXG4gICAqXG4gICAqIEBwYXJhbSB7bnVtYmVyfE9iamVjdHxzdHJpbmd9IHcgcmVhbFxuICAgKiBAcGFyYW0ge251bWJlcj19IHggaW1hZ1xuICAgKiBAcGFyYW0ge251bWJlcj19IHkgaW1hZ1xuICAgKiBAcGFyYW0ge251bWJlcj19IHogaW1hZ1xuICAgKiBAcmV0dXJucyB7UXVhdGVybmlvbn1cbiAgICovXG4gICdwb3cnOiBmdW5jdGlvbiAodywgeCwgeSwgeikge1xuXG4gICAgcGFyc2UoUCwgdywgeCwgeSwgeik7XG5cbiAgICBpZiAoUFsneSddID09PSAwICYmIFBbJ3onXSA9PT0gMCkge1xuXG4gICAgICBpZiAoUFsndyddID09PSAxICYmIFBbJ3gnXSA9PT0gMCkge1xuICAgICAgICByZXR1cm4gdGhpcztcbiAgICAgIH1cblxuICAgICAgaWYgKFBbJ3cnXSA9PT0gMCAmJiBQWyd4J10gPT09IDApIHtcbiAgICAgICAgcmV0dXJuIFF1YXRlcm5pb25bJ09ORSddO1xuICAgICAgfVxuXG4gICAgICAvLyBDaGVjayBpZiB3ZSBjYW4gb3BlcmF0ZSBpbiBDXG4gICAgICAvLyBCb3Jyb3dlZCBmcm9tIGNvbXBsZXguanNcbiAgICAgIGlmICh0aGlzWyd5J10gPT09IDAgJiYgdGhpc1sneiddID09PSAwKSB7XG5cbiAgICAgICAgbGV0IGEgPSB0aGlzWyd3J107XG4gICAgICAgIGxldCBiID0gdGhpc1sneCddO1xuXG4gICAgICAgIGlmIChhID09PSAwICYmIGIgPT09IDApIHtcbiAgICAgICAgICByZXR1cm4gUXVhdGVybmlvblsnWkVSTyddO1xuICAgICAgICB9XG5cbiAgICAgICAgbGV0IGFyZyA9IE1hdGguYXRhbjIoYiwgYSk7XG4gICAgICAgIGxldCBsb2ggPSBsb2dIeXBvdChhLCBiKTtcblxuICAgICAgICBpZiAoUFsneCddID09PSAwKSB7XG5cbiAgICAgICAgICBpZiAoYiA9PT0gMCAmJiBhID49IDApIHtcblxuICAgICAgICAgICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oTWF0aC5wb3coYSwgUFsndyddKSwgMCwgMCwgMCk7XG5cbiAgICAgICAgICB9IGVsc2UgaWYgKGEgPT09IDApIHtcblxuICAgICAgICAgICAgc3dpdGNoIChQWyd3J10gJSA0KSB7XG4gICAgICAgICAgICAgIGNhc2UgMDpcbiAgICAgICAgICAgICAgICByZXR1cm4gbmV3UXVhdGVybmlvbihNYXRoLnBvdyhiLCBQWyd3J10pLCAwLCAwLCAwKTtcbiAgICAgICAgICAgICAgY2FzZSAxOlxuICAgICAgICAgICAgICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKDAsIE1hdGgucG93KGIsIFBbJ3cnXSksIDAsIDApO1xuICAgICAgICAgICAgICBjYXNlIDI6XG4gICAgICAgICAgICAgICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oLU1hdGgucG93KGIsIFBbJ3cnXSksIDAsIDAsIDApO1xuICAgICAgICAgICAgICBjYXNlIDM6XG4gICAgICAgICAgICAgICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oMCwgLU1hdGgucG93KGIsIFBbJ3cnXSksIDAsIDApO1xuICAgICAgICAgICAgfVxuICAgICAgICAgIH1cbiAgICAgICAgfVxuXG4gICAgICAgIGEgPSBNYXRoLmV4cChQWyd3J10gKiBsb2ggLSBQWyd4J10gKiBhcmcpO1xuICAgICAgICBiID0gUFsneCddICogbG9oICsgUFsndyddICogYXJnO1xuICAgICAgICByZXR1cm4gbmV3UXVhdGVybmlvbihcbiAgICAgICAgICBhICogTWF0aC5jb3MoYiksXG4gICAgICAgICAgYSAqIE1hdGguc2luKGIpLCAwLCAwKTtcbiAgICAgIH1cbiAgICB9XG5cbiAgICAvLyBOb3JtYWwgcXVhdGVybmlvbiBiZWhhdmlvclxuICAgIC8vIHFecCA9IGVebG4ocV5wKSA9IGVeKGxuKHEpKnApXG4gICAgcmV0dXJuIHRoaXNbJ2xvZyddKClbJ211bCddKFApWydleHAnXSgpO1xuICB9LFxuICAvKipcbiAgICogQ2hlY2tzIGlmIHR3byBxdWF0cyBhcmUgdGhlIHNhbWVcbiAgICpcbiAgICogQHBhcmFtIHtudW1iZXJ8T2JqZWN0fHN0cmluZ30gdyByZWFsXG4gICAqIEBwYXJhbSB7bnVtYmVyPX0geCBpbWFnXG4gICAqIEBwYXJhbSB7bnVtYmVyPX0geSBpbWFnXG4gICAqIEBwYXJhbSB7bnVtYmVyPX0geiBpbWFnXG4gICAqIEByZXR1cm5zIHtib29sZWFufVxuICAgKi9cbiAgJ2VxdWFscyc6IGZ1bmN0aW9uICh3LCB4LCB5LCB6KSB7XG5cbiAgICBwYXJzZShQLCB3LCB4LCB5LCB6KTtcblxuICAgIGNvbnN0IGVwcyA9IEVQU0lMT047XG5cbiAgICAvLyBtYXliZSBjaGVjayBmb3IgTmFOJ3MgaGVyZT9cbiAgICByZXR1cm4gTWF0aC5hYnMoUFsndyddIC0gdGhpc1sndyddKSA8IGVwc1xuICAgICAgJiYgTWF0aC5hYnMoUFsneCddIC0gdGhpc1sneCddKSA8IGVwc1xuICAgICAgJiYgTWF0aC5hYnMoUFsneSddIC0gdGhpc1sneSddKSA8IGVwc1xuICAgICAgJiYgTWF0aC5hYnMoUFsneiddIC0gdGhpc1sneiddKSA8IGVwcztcbiAgfSxcbiAgLyoqXG4gICAqIENoZWNrcyBpZiBhbGwgcGFydHMgb2YgYSBxdWF0ZXJuaW9uIGFyZSBmaW5pdGVcbiAgICpcbiAgICogQHJldHVybnMge2Jvb2xlYW59XG4gICAqL1xuICAnaXNGaW5pdGUnOiBmdW5jdGlvbiAoKSB7XG5cbiAgICByZXR1cm4gaXNGaW5pdGUodGhpc1sndyddKSAmJiBpc0Zpbml0ZSh0aGlzWyd4J10pICYmIGlzRmluaXRlKHRoaXNbJ3knXSkgJiYgaXNGaW5pdGUodGhpc1sneiddKTtcbiAgfSxcbiAgLyoqXG4gICAqIENoZWNrcyBpZiBhbnkgb2YgdGhlIHBhcnRzIG9mIHRoZSBxdWF0ZXJuaW9uIGlzIG5vdCBhIG51bWJlclxuICAgKlxuICAgKiBAcmV0dXJucyB7Ym9vbGVhbn1cbiAgICovXG4gICdpc05hTic6IGZ1bmN0aW9uICgpIHtcblxuICAgIHJldHVybiBpc05hTih0aGlzWyd3J10pIHx8IGlzTmFOKHRoaXNbJ3gnXSkgfHwgaXNOYU4odGhpc1sneSddKSB8fCBpc05hTih0aGlzWyd6J10pO1xuICB9LFxuICAvKipcbiAgICogR2V0cyB0aGUgUXVhdGVybmlvbiBhcyBhIHdlbGwgZm9ybWF0dGVkIHN0cmluZ1xuICAgKlxuICAgKiBAcmV0dXJucyB7c3RyaW5nfVxuICAgKi9cbiAgJ3RvU3RyaW5nJzogZnVuY3Rpb24gKCkge1xuXG4gICAgY29uc3QgdyA9IHRoaXNbJ3cnXTtcbiAgICBjb25zdCB4ID0gdGhpc1sneCddO1xuICAgIGNvbnN0IHkgPSB0aGlzWyd5J107XG4gICAgY29uc3QgeiA9IHRoaXNbJ3onXTtcbiAgICBsZXQgcmV0ID0gJyc7XG5cbiAgICBpZiAoaXNOYU4odykgfHwgaXNOYU4oeCkgfHwgaXNOYU4oeSkgfHwgaXNOYU4oeikpIHtcbiAgICAgIHJldHVybiAnTmFOJztcbiAgICB9XG5cbiAgICAvLyBBbHRlcm5hdGl2ZSBkZXNpZ24/XG4gICAgLy8gJyglZiwgWyVmICVmICVmXSknXG5cbiAgICByZXQgPSBudW1Ub1N0cih3LCAnJywgcmV0KTtcbiAgICByZXQgKz0gbnVtVG9TdHIoeCwgJ2knLCByZXQpO1xuICAgIHJldCArPSBudW1Ub1N0cih5LCAnaicsIHJldCk7XG4gICAgcmV0ICs9IG51bVRvU3RyKHosICdrJywgcmV0KTtcblxuICAgIGlmICgnJyA9PT0gcmV0KVxuICAgICAgcmV0dXJuICcwJztcblxuICAgIHJldHVybiByZXQ7XG4gIH0sXG4gIC8qKlxuICAgKiBSZXR1cm5zIHRoZSByZWFsIHBhcnQgb2YgdGhlIHF1YXRlcm5pb25cbiAgICpcbiAgICogQHJldHVybnMge251bWJlcn1cbiAgICovXG4gICdyZWFsJzogZnVuY3Rpb24gKCkge1xuXG4gICAgcmV0dXJuIHRoaXNbJ3cnXTtcbiAgfSxcbiAgLyoqXG4gICAqIFJldHVybnMgdGhlIGltYWdpbmFyeSBwYXJ0IG9mIHRoZSBxdWF0ZXJuaW9uIGFzIGEgM0QgdmVjdG9yIC8gYXJyYXlcbiAgICpcbiAgICogQHJldHVybnMge0FycmF5fVxuICAgKi9cbiAgJ2ltYWcnOiBmdW5jdGlvbiAoKSB7XG5cbiAgICByZXR1cm4gW3RoaXNbJ3gnXSwgdGhpc1sneSddLCB0aGlzWyd6J11dO1xuICB9LFxuICAvKipcbiAgICogR2V0cyB0aGUgYWN0dWFsIHF1YXRlcm5pb24gYXMgYSA0RCB2ZWN0b3IgLyBhcnJheVxuICAgKlxuICAgKiBAcmV0dXJucyB7QXJyYXl9XG4gICAqL1xuICAndG9WZWN0b3InOiBmdW5jdGlvbiAoKSB7XG5cbiAgICByZXR1cm4gW3RoaXNbJ3cnXSwgdGhpc1sneCddLCB0aGlzWyd5J10sIHRoaXNbJ3onXV07XG4gIH0sXG4gIC8qKlxuICAgKiBDYWxjdWxhdGVzIHRoZSAzeDMgcm90YXRpb24gbWF0cml4IGZvciB0aGUgY3VycmVudCBxdWF0XG4gICAqXG4gICAqIEBwYXJhbSB7Ym9vbGVhbj19IHR3b0RcbiAgICogQHJldHVybnMge0FycmF5fVxuICAgKi9cbiAgJ3RvTWF0cml4JzogZnVuY3Rpb24gKHR3b0QpIHtcblxuICAgIGNvbnN0IHcgPSB0aGlzWyd3J107XG4gICAgY29uc3QgeCA9IHRoaXNbJ3gnXTtcbiAgICBjb25zdCB5ID0gdGhpc1sneSddO1xuICAgIGNvbnN0IHogPSB0aGlzWyd6J107XG5cbiAgICBjb25zdCB3eCA9IHcgKiB4LCB3eSA9IHcgKiB5LCB3eiA9IHcgKiB6O1xuICAgIGNvbnN0IHh4ID0geCAqIHgsIHh5ID0geCAqIHksIHh6ID0geCAqIHo7XG4gICAgY29uc3QgeXkgPSB5ICogeSwgeXogPSB5ICogeiwgenogPSB6ICogejtcblxuICAgIGlmICh0d29EKSB7XG4gICAgICByZXR1cm4gW1xuICAgICAgICBbMSAtIDIgKiAoeXkgKyB6eiksIDIgKiAoeHkgLSB3eiksIDIgKiAoeHogKyB3eSldLFxuICAgICAgICBbMiAqICh4eSArIHd6KSwgMSAtIDIgKiAoeHggKyB6eiksIDIgKiAoeXogLSB3eCldLFxuICAgICAgICBbMiAqICh4eiAtIHd5KSwgMiAqICh5eiArIHd4KSwgMSAtIDIgKiAoeHggKyB5eSldXTtcbiAgICB9XG5cbiAgICByZXR1cm4gW1xuICAgICAgMSAtIDIgKiAoeXkgKyB6eiksIDIgKiAoeHkgLSB3eiksIDIgKiAoeHogKyB3eSksXG4gICAgICAyICogKHh5ICsgd3opLCAxIC0gMiAqICh4eCArIHp6KSwgMiAqICh5eiAtIHd4KSxcbiAgICAgIDIgKiAoeHogLSB3eSksIDIgKiAoeXogKyB3eCksIDEgLSAyICogKHh4ICsgeXkpXTtcbiAgfSxcbiAgLyoqXG4gICAqIENhbGN1bGF0ZXMgdGhlIGhvbW9nZW5lb3VzIDR4NCByb3RhdGlvbiBtYXRyaXggZm9yIHRoZSBjdXJyZW50IHF1YXRcbiAgICpcbiAgICogQHBhcmFtIHtib29sZWFuPX0gdHdvRFxuICAgKiBAcmV0dXJucyB7QXJyYXl9XG4gICAqL1xuICAndG9NYXRyaXg0JzogZnVuY3Rpb24gKHR3b0QpIHtcblxuICAgIGNvbnN0IHcgPSB0aGlzWyd3J107XG4gICAgY29uc3QgeCA9IHRoaXNbJ3gnXTtcbiAgICBjb25zdCB5ID0gdGhpc1sneSddO1xuICAgIGNvbnN0IHogPSB0aGlzWyd6J107XG5cbiAgICBjb25zdCB3eCA9IHcgKiB4LCB3eSA9IHcgKiB5LCB3eiA9IHcgKiB6O1xuICAgIGNvbnN0IHh4ID0geCAqIHgsIHh5ID0geCAqIHksIHh6ID0geCAqIHo7XG4gICAgY29uc3QgeXkgPSB5ICogeSwgeXogPSB5ICogeiwgenogPSB6ICogejtcblxuICAgIGlmICh0d29EKSB7XG4gICAgICByZXR1cm4gW1xuICAgICAgICBbMSAtIDIgKiAoeXkgKyB6eiksIDIgKiAoeHkgLSB3eiksIDIgKiAoeHogKyB3eSksIDBdLFxuICAgICAgICBbMiAqICh4eSArIHd6KSwgMSAtIDIgKiAoeHggKyB6eiksIDIgKiAoeXogLSB3eCksIDBdLFxuICAgICAgICBbMiAqICh4eiAtIHd5KSwgMiAqICh5eiArIHd4KSwgMSAtIDIgKiAoeHggKyB5eSksIDBdLFxuICAgICAgICBbMCwgMCwgMCwgMV1dO1xuICAgIH1cblxuICAgIHJldHVybiBbXG4gICAgICAxIC0gMiAqICh5eSArIHp6KSwgMiAqICh4eSAtIHd6KSwgMiAqICh4eiArIHd5KSwgMCxcbiAgICAgIDIgKiAoeHkgKyB3eiksIDEgLSAyICogKHh4ICsgenopLCAyICogKHl6IC0gd3gpLCAwLFxuICAgICAgMiAqICh4eiAtIHd5KSwgMiAqICh5eiArIHd4KSwgMSAtIDIgKiAoeHggKyB5eSksIDAsXG4gICAgICAwLCAwLCAwLCAxXTtcbiAgfSxcbiAgLyoqXG4gICAqIERldGVybWluZXMgdGhlIGhvbW9nZW5lb3VzIHJvdGF0aW9uIG1hdHJpeCBzdHJpbmcgdXNlZCBmb3IgQ1NTIDNEIHRyYW5zZm9ybXNcbiAgICpcbiAgICogQHJldHVybnMge3N0cmluZ31cbiAgICovXG4gICd0b0NTU1RyYW5zZm9ybSc6IGZ1bmN0aW9uICgpIHtcblxuICAgIGNvbnN0IHcgPSB0aGlzWyd3J107XG5cbiAgICBsZXQgYW5nbGUgPSAyICogTWF0aC5hY29zKHcpO1xuICAgIGxldCBzaW4yID0gMSAtIHcgKiB3O1xuXG4gICAgaWYgKHNpbjIgPCBFUFNJTE9OKSB7XG4gICAgICBhbmdsZSA9IDA7XG4gICAgICBzaW4yID0gMTtcbiAgICB9IGVsc2Uge1xuICAgICAgc2luMiA9IDEgLyBNYXRoLnNxcnQoc2luMik7IC8vIFJlLXVzZSB2YXJpYWJsZSBzaW5eMiBmb3IgMSAvIHNpblxuICAgIH1cbiAgICByZXR1cm4gXCJyb3RhdGUzZChcIiArIHRoaXNbJ3gnXSAqIHNpbjIgKyBcIixcIiArIHRoaXNbJ3knXSAqIHNpbjIgKyBcIixcIiArIHRoaXNbJ3onXSAqIHNpbjIgKyBcIixcIiArIGFuZ2xlICsgXCJyYWQpXCI7XG4gIH0sXG4gIC8qKlxuICAgKiBDYWxjdWxhdGVzIHRoZSBheGlzIGFuZCBhbmdsZSByZXByZXNlbnRhdGlvbiBvZiB0aGUgcXVhdGVybmlvblxuICAgKlxuICAgKiBAcmV0dXJucyB7QXJyYXl9XG4gICAqL1xuICAndG9BeGlzQW5nbGUnOiBmdW5jdGlvbiAoKSB7XG5cbiAgICBjb25zdCB3ID0gdGhpc1sndyddO1xuICAgIGNvbnN0IHNpbjIgPSAxIC0gdyAqIHc7IC8vIHNpbihhbmdsZSAvIDIpID0gc2luKGFjb3ModykpID0gc3FydCgxIC0gd14yKSA9IHx2fCwgc2luY2UgMSA9IGRvdChRKSA8PT4gZG90KHYpID0gMSAtIHdeMlxuXG4gICAgaWYgKHNpbjIgPCBFUFNJTE9OKSB7IC8vIEFsdGVybmF0aXZlbHkgfHZ8ID09IDBcbiAgICAgIC8vIElmIHRoZSBzaW5lIGlzIGNsb3NlIHRvIDAsIHdlJ3JlIGNsb3NlIHRvIHRoZSB1bml0IHF1YXRlcm5pb24gYW5kIHRoZSBkaXJlY3Rpb24gZG9lcyBub3QgbWF0dGVyXG4gICAgICByZXR1cm4gW1t0aGlzWyd4J10sIHRoaXNbJ3knXSwgdGhpc1sneiddXSwgMF07IC8vIG9yIFtbMSwgMCwgMF0sIDBdID8gIG9yIFtbMCwgMCwgMF0sIDBdID9cbiAgICB9XG5cbiAgICBjb25zdCBpc2luID0gMSAvIE1hdGguc3FydChzaW4yKTtcbiAgICBjb25zdCBhbmdsZSA9IDIgKiBNYXRoLmFjb3Modyk7IC8vIEFsdGVybmF0aXZlbHk6IDIgKiBhdGFuMih8dnwsIHcpXG4gICAgcmV0dXJuIFtbdGhpc1sneCddICogaXNpbiwgdGhpc1sneSddICogaXNpbiwgdGhpc1sneiddICogaXNpbl0sIGFuZ2xlXTtcbiAgfSxcbiAgLyoqXG4gICAqIENhbGN1bGF0ZXMgdGhlIEV1bGVyIGFuZ2xlcyByZXByZXNlbnRlZCBieSB0aGUgY3VycmVudCBxdWF0IChtdWx0aXBsaWNhdGlvbiBvcmRlciBmcm9tIHJpZ2h0IHRvIGxlZnQpXG4gICAqIFxuICAgKiBAcGFyYW0ge3N0cmluZz19IG9yZGVyIEF4aXMgb3JkZXIgKFRhaXQgQnJ5YW4pXG4gICAqIEByZXR1cm5zIHtPYmplY3R9XG4gICAqL1xuICAndG9FdWxlcic6IGZ1bmN0aW9uIChvcmRlcikge1xuXG4gICAgY29uc3QgdyA9IHRoaXNbJ3cnXTtcbiAgICBjb25zdCB4ID0gdGhpc1sneCddO1xuICAgIGNvbnN0IHkgPSB0aGlzWyd5J107XG4gICAgY29uc3QgeiA9IHRoaXNbJ3onXTtcblxuICAgIGNvbnN0IHd4ID0gdyAqIHgsIHd5ID0gdyAqIHksIHd6ID0gdyAqIHo7XG4gICAgY29uc3QgeHggPSB4ICogeCwgeHkgPSB4ICogeSwgeHogPSB4ICogejtcbiAgICBjb25zdCB5eSA9IHkgKiB5LCB5eiA9IHkgKiB6LCB6eiA9IHogKiB6O1xuXG4gICAgZnVuY3Rpb24gYXNpbih0KSB7XG4gICAgICByZXR1cm4gdCA+PSAxID8gTWF0aC5QSSAvIDIgOiAodCA8PSAtMSA/IC1NYXRoLlBJIC8gMiA6IE1hdGguYXNpbih0KSk7XG4gICAgfVxuXG4gICAgaWYgKG9yZGVyID09PSB1bmRlZmluZWQgfHwgb3JkZXIgPT09ICdaWFknKSB7XG4gICAgICByZXR1cm4gW1xuICAgICAgICAtTWF0aC5hdGFuMigyICogKHh5IC0gd3opLCAxIC0gMiAqICh4eCArIHp6KSksXG4gICAgICAgIGFzaW4oMiAqICh5eiArIHd4KSksXG4gICAgICAgIC1NYXRoLmF0YW4yKDIgKiAoeHogLSB3eSksIDEgLSAyICogKHh4ICsgeXkpKSxcbiAgICAgIF07XG4gICAgfVxuXG4gICAgaWYgKG9yZGVyID09PSAnWFlaJyB8fCBvcmRlciA9PT0gJ1JQWScpIHtcbiAgICAgIHJldHVybiBbXG4gICAgICAgIC1NYXRoLmF0YW4yKDIgKiAoeXogLSB3eCksIDEgLSAyICogKHh4ICsgeXkpKSxcbiAgICAgICAgYXNpbigyICogKHh6ICsgd3kpKSxcbiAgICAgICAgLU1hdGguYXRhbjIoMiAqICh4eSAtIHd6KSwgMSAtIDIgKiAoeXkgKyB6eikpLFxuICAgICAgXTtcbiAgICB9XG5cbiAgICBpZiAob3JkZXIgPT09ICdZWFonKSB7XG4gICAgICByZXR1cm4gW1xuICAgICAgICBNYXRoLmF0YW4yKDIgKiAoeHogKyB3eSksIDEgLSAyICogKHh4ICsgeXkpKSxcbiAgICAgICAgLWFzaW4oMiAqICh5eiAtIHd4KSksXG4gICAgICAgIE1hdGguYXRhbjIoMiAqICh4eSArIHd6KSwgMSAtIDIgKiAoeHggKyB6eikpLFxuICAgICAgXTtcbiAgICB9XG5cbiAgICBpZiAob3JkZXIgPT09ICdaWVgnIHx8IG9yZGVyID09PSAnWVBSJykgeyAgLy8gcm9sbCBhcm91bmQgWCwgcGl0Y2ggYXJvdW5kIFksIHlhdyBhcm91bmQgWlxuICAgICAgLypcbiAgICAgIGlmICgyICogKHh6IC0gd3kpID4gLjk5OSkge1xuICAgICAgICByZXR1cm4gW1xuICAgICAgICAgIDIgKiBNYXRoLmF0YW4yKHgsIHcpLFxuICAgICAgICAgIC1NYXRoLlBJIC8gMixcbiAgICAgICAgICAwXG4gICAgICAgIF07XG4gICAgICB9XG5cbiAgICAgIGlmICgyICogKHh6IC0gd3kpIDwgLS45OTkpIHtcbiAgICAgICAgcmV0dXJuIFtcbiAgICAgICAgICAtMiAqIE1hdGguYXRhbjIoeCwgdyksXG4gICAgICAgICAgTWF0aC5QSSAvIDIsXG4gICAgICAgICAgMFxuICAgICAgICBdO1xuICAgICAgfVxuICAgICAgKi9cbiAgICAgIHJldHVybiBbXG4gICAgICAgIE1hdGguYXRhbjIoMiAqICh4eSArIHd6KSwgMSAtIDIgKiAoeXkgKyB6eikpLCAvLyBIZWFkaW5nIC8gWWF3XG4gICAgICAgIC1hc2luKDIgKiAoeHogLSB3eSkpLCAvLyBBdHRpdHVkZSAvIFBpdGNoXG4gICAgICAgIE1hdGguYXRhbjIoMiAqICh5eiArIHd4KSwgMSAtIDIgKiAoeHggKyB5eSkpLCAvLyBCYW5rIC8gUm9sbFxuICAgICAgXTtcbiAgICB9XG5cbiAgICBpZiAob3JkZXIgPT09ICdZWlgnKSB7XG4gICAgICAvKlxuICAgICAgaWYgKDIgKiAoeHkgKyB3eikgPiAuOTk5KSB7IC8vIE5vcnRoIHBvbGVcbiAgICAgICAgcmV0dXJuIFtcbiAgICAgICAgICAyICogTWF0aC5hdGFuMih4LCB3KSxcbiAgICAgICAgICBNYXRoLlBJIC8gMixcbiAgICAgICAgICAwXG4gICAgICAgIF07XG4gICAgICB9XG5cbiAgICAgIGlmICgyICogKHh5ICsgd3opIDwgLS45OTkpIHsgLy8gU291dGggcG9sZVxuICAgICAgICByZXR1cm4gW1xuICAgICAgICAgIC0yICogTWF0aC5hdGFuMih4LCB3KSxcbiAgICAgICAgICAtTWF0aC5QSSAvIDIsXG4gICAgICAgICAgMFxuICAgICAgICBdO1xuICAgICAgfVxuICAgICAgKi9cbiAgICAgIHJldHVybiBbXG4gICAgICAgIC1NYXRoLmF0YW4yKDIgKiAoeHogLSB3eSksIDEgLSAyICogKHl5ICsgenopKSwgLy8gSGVhZGluZ1xuICAgICAgICBhc2luKDIgKiAoeHkgKyB3eikpLCAvLyBBdHRpdHVkZVxuICAgICAgICAtTWF0aC5hdGFuMigyICogKHl6IC0gd3gpLCAxIC0gMiAqICh4eCArIHp6KSksIC8vIEJhbmtcbiAgICAgIF07XG4gICAgfVxuXG4gICAgaWYgKG9yZGVyID09PSAnWFpZJykge1xuICAgICAgcmV0dXJuIFtcbiAgICAgICAgTWF0aC5hdGFuMigyICogKHl6ICsgd3gpLCAxIC0gMiAqICh4eCArIHp6KSksXG4gICAgICAgIC1hc2luKDIgKiAoeHkgLSB3eikpLFxuICAgICAgICBNYXRoLmF0YW4yKDIgKiAoeHogKyB3eSksIDEgLSAyICogKHl5ICsgenopKSxcbiAgICAgIF07XG4gICAgfVxuICAgIHJldHVybiBudWxsO1xuICB9LFxuICAvKipcbiAgICogQ2xvbmVzIHRoZSBhY3R1YWwgb2JqZWN0XG4gICAqXG4gICAqIEByZXR1cm5zIHtRdWF0ZXJuaW9ufVxuICAgKi9cbiAgJ2Nsb25lJzogZnVuY3Rpb24gKCkge1xuXG4gICAgcmV0dXJuIG5ld1F1YXRlcm5pb24odGhpc1sndyddLCB0aGlzWyd4J10sIHRoaXNbJ3knXSwgdGhpc1sneiddKTtcbiAgfSxcbiAgLyoqXG4gICAqIFJvdGF0ZXMgYSB2ZWN0b3IgYWNjb3JkaW5nIHRvIHRoZSBjdXJyZW50IHF1YXRlcm5pb24sIGFzc3VtZXMgfHF8PTFcbiAgICogQGxpbmsgaHR0cHM6Ly9yYXcub3JnL3Byb29mL3ZlY3Rvci1yb3RhdGlvbi11c2luZy1xdWF0ZXJuaW9ucy9cbiAgICpcbiAgICogQHBhcmFtIHtBcnJheX0gdiBUaGUgdmVjdG9yIHRvIGJlIHJvdGF0ZWRcbiAgICogQHJldHVybnMge0FycmF5fVxuICAgKi9cbiAgJ3JvdGF0ZVZlY3Rvcic6IGZ1bmN0aW9uICh2KSB7XG5cbiAgICBjb25zdCBxdyA9IHRoaXNbJ3cnXTtcbiAgICBjb25zdCBxeCA9IHRoaXNbJ3gnXTtcbiAgICBjb25zdCBxeSA9IHRoaXNbJ3knXTtcbiAgICBjb25zdCBxeiA9IHRoaXNbJ3onXTtcblxuICAgIGNvbnN0IHZ4ID0gdlswXTtcbiAgICBjb25zdCB2eSA9IHZbMV07XG4gICAgY29uc3QgdnogPSB2WzJdO1xuXG4gICAgLy8gdCA9IHEgeCB2XG4gICAgbGV0IHR4ID0gcXkgKiB2eiAtIHF6ICogdnk7XG4gICAgbGV0IHR5ID0gcXogKiB2eCAtIHF4ICogdno7XG4gICAgbGV0IHR6ID0gcXggKiB2eSAtIHF5ICogdng7XG5cbiAgICAvLyB0ID0gMnRcbiAgICB0eCA9IHR4ICsgdHg7XG4gICAgdHkgPSB0eSArIHR5O1xuICAgIHR6ID0gdHogKyB0ejtcblxuICAgIC8vIHYgKyB3IHQgKyBxIHggdFxuICAgIHJldHVybiBbXG4gICAgICB2eCArIHF3ICogdHggKyBxeSAqIHR6IC0gcXogKiB0eSxcbiAgICAgIHZ5ICsgcXcgKiB0eSArIHF6ICogdHggLSBxeCAqIHR6LFxuICAgICAgdnogKyBxdyAqIHR6ICsgcXggKiB0eSAtIHF5ICogdHhdO1xuICB9LFxuXG4gIC8qKlxuICAgKiBHZXRzIGEgZnVuY3Rpb24gdG8gc3BoZXJpY2FsbHkgaW50ZXJwb2xhdGUgYmV0d2VlbiB0d28gcXVhdGVybmlvbnNcbiAgICogXG4gICAqIEByZXR1cm5zIEZ1bmN0aW9uXG4gICAqL1xuICAnc2xlcnAnOiBmdW5jdGlvbiAodywgeCwgeSwgeikge1xuXG4gICAgcGFyc2UoUCwgdywgeCwgeSwgeik7XG5cbiAgICAvLyBzbGVycChRMSwgUTIsIHQpIDo9IFExKFExXi0xIFEyKV50XG5cbiAgICBsZXQgdzEgPSB0aGlzWyd3J107XG4gICAgbGV0IHgxID0gdGhpc1sneCddO1xuICAgIGxldCB5MSA9IHRoaXNbJ3knXTtcbiAgICBsZXQgejEgPSB0aGlzWyd6J107XG5cbiAgICBsZXQgdzIgPSBQWyd3J107XG4gICAgbGV0IHgyID0gUFsneCddO1xuICAgIGxldCB5MiA9IFBbJ3knXTtcbiAgICBsZXQgejIgPSBQWyd6J107XG5cbiAgICBsZXQgY29zVGhldGEwID0gdzEgKiB3MiArIHgxICogeDIgKyB5MSAqIHkyICsgejEgKiB6MjtcblxuICAgIGlmIChjb3NUaGV0YTAgPCAwKSB7XG4gICAgICB3MSA9IC13MTtcbiAgICAgIHgxID0gLXgxO1xuICAgICAgeTEgPSAteTE7XG4gICAgICB6MSA9IC16MTtcbiAgICAgIGNvc1RoZXRhMCA9IC1jb3NUaGV0YTA7XG4gICAgfVxuXG4gICAgaWYgKGNvc1RoZXRhMCA+PSAxIC0gRVBTSUxPTikge1xuICAgICAgcmV0dXJuIGZ1bmN0aW9uIChwY3QpIHtcbiAgICAgICAgcmV0dXJuIG5ld05vcm1hbGl6ZWQoXG4gICAgICAgICAgdzEgKyBwY3QgKiAodzIgLSB3MSksXG4gICAgICAgICAgeDEgKyBwY3QgKiAoeDIgLSB4MSksXG4gICAgICAgICAgeTEgKyBwY3QgKiAoeTIgLSB5MSksXG4gICAgICAgICAgejEgKyBwY3QgKiAoejIgLSB6MSkpO1xuICAgICAgfTtcbiAgICB9XG5cbiAgICBsZXQgVGhldGEwID0gTWF0aC5hY29zKGNvc1RoZXRhMCk7XG4gICAgbGV0IHNpblRoZXRhMCA9IE1hdGguc2luKFRoZXRhMCk7XG5cbiAgICByZXR1cm4gZnVuY3Rpb24gKHBjdCkge1xuXG4gICAgICBsZXQgVGhldGEgPSBUaGV0YTAgKiBwY3Q7XG4gICAgICBsZXQgc2luVGhldGEgPSBNYXRoLnNpbihUaGV0YSk7XG4gICAgICBsZXQgY29zVGhldGEgPSBNYXRoLmNvcyhUaGV0YSk7XG5cbiAgICAgIGxldCBzMCA9IGNvc1RoZXRhIC0gY29zVGhldGEwICogc2luVGhldGEgLyBzaW5UaGV0YTA7XG4gICAgICBsZXQgczEgPSBzaW5UaGV0YSAvIHNpblRoZXRhMDtcblxuICAgICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oXG4gICAgICAgIHMwICogdzEgKyBzMSAqIHcyLFxuICAgICAgICBzMCAqIHgxICsgczEgKiB4MixcbiAgICAgICAgczAgKiB5MSArIHMxICogeTIsXG4gICAgICAgIHMwICogejEgKyBzMSAqIHoyKTtcbiAgICB9O1xuICB9XG59O1xuXG5RdWF0ZXJuaW9uWydaRVJPJ10gPSBuZXdRdWF0ZXJuaW9uKDAsIDAsIDAsIDApOyAvLyBUaGlzIGlzIHRoZSBhZGRpdGl2ZSBpZGVudGl0eSBRdWF0ZXJuaW9uXG5RdWF0ZXJuaW9uWydPTkUnXSA9IG5ld1F1YXRlcm5pb24oMSwgMCwgMCwgMCk7IC8vIFRoaXMgaXMgdGhlIG11bHRpcGxpY2F0aXZlIGlkZW50aXR5IFF1YXRlcm5pb25cblF1YXRlcm5pb25bJ0knXSA9IG5ld1F1YXRlcm5pb24oMCwgMSwgMCwgMCk7XG5RdWF0ZXJuaW9uWydKJ10gPSBuZXdRdWF0ZXJuaW9uKDAsIDAsIDEsIDApO1xuUXVhdGVybmlvblsnSyddID0gbmV3UXVhdGVybmlvbigwLCAwLCAwLCAxKTtcblxuLyoqXG4gKiBAY29uc3RcbiAqL1xuY29uc3QgRVBTSUxPTiA9IDFlLTE2O1xuXG4vKipcbiAqIENyZWF0ZXMgcXVhdGVybmlvbiBieSBhIHJvdGF0aW9uIGdpdmVuIGFzIGF4aXMtYW5nbGUgb3JpZW50YXRpb25cbiAqXG4gKiBAcGFyYW0ge0FycmF5fSBheGlzIFRoZSBheGlzIGFyb3VuZCB3aGljaCB0byByb3RhdGVcbiAqIEBwYXJhbSB7bnVtYmVyfSBhbmdsZSBUaGUgYW5nbGUgaW4gcmFkaWFuc1xuICogQHJldHVybnMge1F1YXRlcm5pb259XG4gKi9cblF1YXRlcm5pb25bJ2Zyb21BeGlzQW5nbGUnXSA9IGZ1bmN0aW9uIChheGlzLCBhbmdsZSkge1xuXG4gIC8vIFEgPSBbY29zKGFuZ2xlIC8gMiksIHYgKiBzaW4oYW5nbGUgLyAyKV1cblxuICBjb25zdCBhID0gYXhpc1swXTtcbiAgY29uc3QgYiA9IGF4aXNbMV07XG4gIGNvbnN0IGMgPSBheGlzWzJdO1xuXG4gIGNvbnN0IGhhbGZBbmdsZSA9IGFuZ2xlICogMC41O1xuXG4gIGNvbnN0IHNpbl8yID0gTWF0aC5zaW4oaGFsZkFuZ2xlKTtcbiAgY29uc3QgY29zXzIgPSBNYXRoLmNvcyhoYWxmQW5nbGUpO1xuXG4gIGNvbnN0IHNpbl9ub3JtID0gc2luXzIgLyBNYXRoLnNxcnQoYSAqIGEgKyBiICogYiArIGMgKiBjKTtcblxuICByZXR1cm4gbmV3UXVhdGVybmlvbihjb3NfMiwgYSAqIHNpbl9ub3JtLCBiICogc2luX25vcm0sIGMgKiBzaW5fbm9ybSk7XG59O1xuXG4vKipcbiAqIENhbGN1bGF0ZXMgdGhlIHF1YXRlcm5pb24gdG8gcm90YXRlIHZlY3RvciB1IG9udG8gdmVjdG9yIHZcbiAqIEBsaW5rIGh0dHBzOi8vcmF3Lm9yZy9wcm9vZi9xdWF0ZXJuaW9uLWZyb20tdHdvLXZlY3RvcnMvXG4gKlxuICogQHBhcmFtIHtBcnJheX0gdVxuICogQHBhcmFtIHtBcnJheX0gdlxuICovXG5RdWF0ZXJuaW9uWydmcm9tVmVjdG9ycyddID0gZnVuY3Rpb24gKHUsIHYpIHtcblxuICBsZXQgdXggPSB1WzBdO1xuICBsZXQgdXkgPSB1WzFdO1xuICBsZXQgdXogPSB1WzJdO1xuXG4gIGxldCB2eCA9IHZbMF07XG4gIGxldCB2eSA9IHZbMV07XG4gIGxldCB2eiA9IHZbMl07XG5cbiAgY29uc3QgdUxlbiA9IE1hdGguc3FydCh1eCAqIHV4ICsgdXkgKiB1eSArIHV6ICogdXopO1xuICBjb25zdCB2TGVuID0gTWF0aC5zcXJ0KHZ4ICogdnggKyB2eSAqIHZ5ICsgdnogKiB2eik7XG5cbiAgLy8gTm9ybWFsaXplIHUgYW5kIHZcbiAgaWYgKHVMZW4gPiAwKSB1eCAvPSB1TGVuLCB1eSAvPSB1TGVuLCB1eiAvPSB1TGVuO1xuICBpZiAodkxlbiA+IDApIHZ4IC89IHZMZW4sIHZ5IC89IHZMZW4sIHZ6IC89IHZMZW47XG5cbiAgLy8gQ2FsY3VsYXRlIGRvdCBwcm9kdWN0IG9mIG5vcm1hbGl6ZWQgdSBhbmQgdlxuICBjb25zdCBkb3QgPSB1eCAqIHZ4ICsgdXkgKiB2eSArIHV6ICogdno7XG5cbiAgLy8gUGFyYWxsZWwgd2hlbiBkb3QgPiAwLjk5OTk5OVxuICBpZiAoZG90ID49IDEgLSBFUFNJTE9OKSB7XG4gICAgcmV0dXJuIFF1YXRlcm5pb25bJ09ORSddO1xuICB9XG5cbiAgLy8gQW50aS1QYXJhbGxlbCAoY2xvc2UgdG8gUEkpIHdoZW4gZG90IDwgLTAuOTk5OTk5XG4gIGlmICgxICsgZG90IDw9IEVQU0lMT04pIHtcblxuICAgIC8vIFJvdGF0ZSAxODDCsCBhcm91bmQgYW55IG9ydGhvZ29uYWwgdmVjdG9yXG4gICAgLy8gYXhpcyA9IGxlbihjcm9zcyhbMSwgMCwgMF0sIHUpKSA9PSAwID8gY3Jvc3MoWzAsIDEsIDBdLCB1KSA6IGNyb3NzKFsxLCAwLCAwXSwgdSkgYW5kIHRoZXJlZm9yZVxuICAgIC8vICAgIHJldHVybiBRdWF0ZXJuaW9uWydmcm9tQXhpc0FuZ2xlJ10oTWF0aC5hYnModXgpID4gTWF0aC5hYnModXopID8gWy11eSwgdXgsIDBdIDogWzAsIC11eiwgdXldLCBNYXRoLlBJKVxuICAgIC8vIG9yIHJldHVybiBRdWF0ZXJuaW9uWydmcm9tQXhpc0FuZ2xlJ10oTWF0aC5hYnModXgpID4gTWF0aC5hYnModXopID8gWyB1eSwtdXgsIDBdIDogWzAsICB1eiwtdXldLCBNYXRoLlBJKVxuICAgIC8vIG9yIC4uLlxuXG4gICAgLy8gU2luY2UgZnJvbUF4aXNBbmdsZShheGlzLCBQSSkgPT0gUXVhdGVybmlvbigwLCBheGlzKS5ub3JtYWxpemUoKSxcbiAgICBpZiAoTWF0aC5hYnModXgpID4gTWF0aC5hYnModXopKSB7XG4gICAgICByZXR1cm4gbmV3Tm9ybWFsaXplZCgwLCAtdXksIHV4LCAwKTtcbiAgICB9IGVsc2Uge1xuICAgICAgcmV0dXJuIG5ld05vcm1hbGl6ZWQoMCwgMCwgLXV6LCB1eSk7XG4gICAgfVxuICB9XG5cbiAgLy8gdyA9IGNyb3NzKHUsIHYpXG4gIGNvbnN0IHd4ID0gdXkgKiB2eiAtIHV6ICogdnk7XG4gIGNvbnN0IHd5ID0gdXogKiB2eCAtIHV4ICogdno7XG4gIGNvbnN0IHd6ID0gdXggKiB2eSAtIHV5ICogdng7XG5cbiAgLy8gfFF8ID0gc3FydCgoMS4wICsgZG90KSAqIDIuMClcbiAgcmV0dXJuIG5ld05vcm1hbGl6ZWQoMSArIGRvdCwgd3gsIHd5LCB3eik7XG59O1xuXG4vKipcbiAqIEdldHMgYSBzcGhlcmljYWwgcmFuZG9tIG51bWJlclxuICogQGxpbmsgaHR0cDovL3BsYW5uaW5nLmNzLnVpdWMuZWR1L25vZGUxOTguaHRtbFxuICovXG5RdWF0ZXJuaW9uWydyYW5kb20nXSA9IGZ1bmN0aW9uICgpIHtcblxuICBjb25zdCB1MSA9IE1hdGgucmFuZG9tKCk7XG4gIGNvbnN0IHUyID0gTWF0aC5yYW5kb20oKTtcbiAgY29uc3QgdTMgPSBNYXRoLnJhbmRvbSgpO1xuXG4gIGNvbnN0IHMgPSBNYXRoLnNxcnQoMSAtIHUxKTtcbiAgY29uc3QgdCA9IE1hdGguc3FydCh1MSk7XG5cbiAgcmV0dXJuIG5ld1F1YXRlcm5pb24oXG4gICAgdCAqIE1hdGguY29zKDIgKiBNYXRoLlBJICogdTMpLFxuICAgIHMgKiBNYXRoLnNpbigyICogTWF0aC5QSSAqIHUyKSxcbiAgICBzICogTWF0aC5jb3MoMiAqIE1hdGguUEkgKiB1MiksXG4gICAgdCAqIE1hdGguc2luKDIgKiBNYXRoLlBJICogdTMpXG4gICk7XG59O1xuXG4vKipcbiAqIENyZWF0ZXMgYSBxdWF0ZXJuaW9uIGJ5IGEgcm90YXRpb24gZ2l2ZW4gYnkgRXVsZXIgYW5nbGVzIChsb2dpY2FsIGFwcGxpY2F0aW9uIG9yZGVyIGZyb20gbGVmdCB0byByaWdodClcbiAqXG4gKiBAcGFyYW0ge251bWJlcn0gz4ggRmlyc3QgYW5nbGVcbiAqIEBwYXJhbSB7bnVtYmVyfSDOuCBTZWNvbmQgYW5nbGVcbiAqIEBwYXJhbSB7bnVtYmVyfSDPhiBUaGlyZCBhbmdsZVxuICogQHBhcmFtIHtzdHJpbmc9fSBvcmRlciBBeGlzIG9yZGVyIChUYWl0IEJyeWFuKVxuICogQHJldHVybnMge1F1YXRlcm5pb259XG4gKi9cblF1YXRlcm5pb25bJ2Zyb21FdWxlckxvZ2ljYWwnXSA9IGZ1bmN0aW9uICjPiCwgzrgsIM+GLCBvcmRlcikge1xuXG4gIHJldHVybiBRdWF0ZXJuaW9uWydmcm9tRXVsZXInXSjPhiwgzrgsIM+ILCBvcmRlciAhPT0gdW5kZWZpbmVkID8gb3JkZXJbMl0gKyBvcmRlclsxXSArIG9yZGVyWzBdIDogb3JkZXIpO1xufTtcblxuLyoqXG4gKiBDcmVhdGVzIGEgcXVhdGVybmlvbiBieSBhIHJvdGF0aW9uIGdpdmVuIGJ5IEV1bGVyIGFuZ2xlcyAobXVsdGlwbGljYXRpb24gb3JkZXIgZnJvbSByaWdodCB0byBsZWZ0KVxuICpcbiAqIEBwYXJhbSB7bnVtYmVyfSDPhiBGaXJzdCBhbmdsZVxuICogQHBhcmFtIHtudW1iZXJ9IM64IFNlY29uZCBhbmdsZVxuICogQHBhcmFtIHtudW1iZXJ9IM+IIFRoaXJkIGFuZ2xlXG4gKiBAcGFyYW0ge3N0cmluZz19IG9yZGVyIEF4aXMgb3JkZXIgKFRhaXQgQnJ5YW4pXG4gKiBAcmV0dXJucyB7UXVhdGVybmlvbn1cbiAqL1xuUXVhdGVybmlvblsnZnJvbUV1bGVyJ10gPSBmdW5jdGlvbiAoz4YsIM64LCDPiCwgb3JkZXIpIHtcblxuICBjb25zdCBfeCA9IM+GICogMC41O1xuICBjb25zdCBfeSA9IM64ICogMC41O1xuICBjb25zdCBfeiA9IM+IICogMC41O1xuXG4gIGNvbnN0IGNYID0gTWF0aC5jb3MoX3gpO1xuICBjb25zdCBjWSA9IE1hdGguY29zKF95KTtcbiAgY29uc3QgY1ogPSBNYXRoLmNvcyhfeik7XG5cbiAgY29uc3Qgc1ggPSBNYXRoLnNpbihfeCk7XG4gIGNvbnN0IHNZID0gTWF0aC5zaW4oX3kpO1xuICBjb25zdCBzWiA9IE1hdGguc2luKF96KTtcblxuICBpZiAob3JkZXIgPT09IHVuZGVmaW5lZCB8fCBvcmRlciA9PT0gJ1pYWScpIHtcbiAgICAvLyBheGlzQW5nbGUoWzAsIDAsIDFdLCDPhikgKiBheGlzQW5nbGUoWzEsIDAsIDBdLCDOuCkgKiBheGlzQW5nbGUoWzAsIDEsIDBdLCDPiClcbiAgICByZXR1cm4gbmV3UXVhdGVybmlvbihcbiAgICAgIGNYICogY1kgKiBjWiAtIHNYICogc1kgKiBzWixcbiAgICAgIHNZICogY1ggKiBjWiAtIHNYICogc1ogKiBjWSxcbiAgICAgIHNYICogc1kgKiBjWiArIHNaICogY1ggKiBjWSxcbiAgICAgIHNYICogY1kgKiBjWiArIHNZICogc1ogKiBjWCk7XG4gIH1cblxuICBpZiAob3JkZXIgPT09ICdYWVonIHx8IG9yZGVyID09PSAnUlBZJykgeyAvLyByb2xsIGFyb3VuZCBYLCBwaXRjaCBhcm91bmQgWSwgeWF3IGFyb3VuZCBaXG4gICAgLy8gYXhpc0FuZ2xlKFsxLCAwLCAwXSwgz4YpICogYXhpc0FuZ2xlKFswLCAxLCAwXSwgzrgpICogYXhpc0FuZ2xlKFswLCAwLCAxXSwgz4gpXG4gICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oXG4gICAgICBjWCAqIGNZICogY1ogLSBzWCAqIHNZICogc1osXG4gICAgICBzWCAqIGNZICogY1ogKyBzWSAqIHNaICogY1gsXG4gICAgICBzWSAqIGNYICogY1ogLSBzWCAqIHNaICogY1ksXG4gICAgICBzWCAqIHNZICogY1ogKyBzWiAqIGNYICogY1kpO1xuICB9XG5cbiAgaWYgKG9yZGVyID09PSAnWVhaJykgeyAvLyBkZXZpY2VvcmllbnRhdGlvblxuICAgIC8vIGF4aXNBbmdsZShbMCwgMSwgMF0sIM+GKSAqIGF4aXNBbmdsZShbMSwgMCwgMF0sIM64KSAqIGF4aXNBbmdsZShbMCwgMCwgMV0sIM+IKVxuICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKFxuICAgICAgc1ggKiBzWSAqIHNaICsgY1ggKiBjWSAqIGNaLFxuICAgICAgc1ggKiBzWiAqIGNZICsgc1kgKiBjWCAqIGNaLFxuICAgICAgc1ggKiBjWSAqIGNaIC0gc1kgKiBzWiAqIGNYLFxuICAgICAgc1ogKiBjWCAqIGNZIC0gc1ggKiBzWSAqIGNaKTtcbiAgfVxuXG4gIGlmIChvcmRlciA9PT0gJ1pZWCcgfHwgb3JkZXIgPT09ICdZUFInKSB7XG4gICAgLy8gYXhpc0FuZ2xlKFswLCAwLCAxXSwgz4YpICogYXhpc0FuZ2xlKFswLCAxLCAwXSwgzrgpICogYXhpc0FuZ2xlKFsxLCAwLCAwXSwgz4gpXG4gICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oXG4gICAgICBzWCAqIHNZICogc1ogKyBjWCAqIGNZICogY1osXG4gICAgICBzWiAqIGNYICogY1kgLSBzWCAqIHNZICogY1osXG4gICAgICBzWCAqIHNaICogY1kgKyBzWSAqIGNYICogY1osXG4gICAgICBzWCAqIGNZICogY1ogLSBzWSAqIHNaICogY1gpO1xuICB9XG5cbiAgaWYgKG9yZGVyID09PSAnWVpYJykge1xuICAgIC8vIGF4aXNBbmdsZShbMCwgMSwgMF0sIM+GKSAqIGF4aXNBbmdsZShbMCwgMCwgMV0sIM64KSAqIGF4aXNBbmdsZShbMSwgMCwgMF0sIM+IKVxuICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKFxuICAgICAgY1ggKiBjWSAqIGNaIC0gc1ggKiBzWSAqIHNaLFxuICAgICAgc1ggKiBzWSAqIGNaICsgc1ogKiBjWCAqIGNZLFxuICAgICAgc1ggKiBjWSAqIGNaICsgc1kgKiBzWiAqIGNYLFxuICAgICAgc1kgKiBjWCAqIGNaIC0gc1ggKiBzWiAqIGNZKTtcbiAgfVxuXG4gIGlmIChvcmRlciA9PT0gJ1haWScpIHtcbiAgICAvLyBheGlzQW5nbGUoWzEsIDAsIDBdLCDPhikgKiBheGlzQW5nbGUoWzAsIDAsIDFdLCDOuCkgKiBheGlzQW5nbGUoWzAsIDEsIDBdLCDPiClcbiAgICByZXR1cm4gbmV3UXVhdGVybmlvbihcbiAgICAgIHNYICogc1kgKiBzWiArIGNYICogY1kgKiBjWixcbiAgICAgIHNYICogY1kgKiBjWiAtIHNZICogc1ogKiBjWCxcbiAgICAgIHNaICogY1ggKiBjWSAtIHNYICogc1kgKiBjWixcbiAgICAgIHNYICogc1ogKiBjWSArIHNZICogY1ggKiBjWik7XG4gIH1cblxuICBpZiAob3JkZXIgPT09ICdaWVonKSB7XG4gICAgLy8gYXhpc0FuZ2xlKFswLCAwLCAxXSwgz4YpICogYXhpc0FuZ2xlKFswLCAxLCAwXSwgzrgpICogYXhpc0FuZ2xlKFswLCAwLCAxXSwgz4gpXG4gICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oXG4gICAgICBjWCAqIGNZICogY1ogLSBzWCAqIHNaICogY1ksXG4gICAgICBzWSAqIHNaICogY1ggLSBzWCAqIHNZICogY1osXG4gICAgICBzWCAqIHNZICogc1ogKyBzWSAqIGNYICogY1osXG4gICAgICBzWCAqIGNZICogY1ogKyBzWiAqIGNYICogY1kpO1xuICB9XG5cbiAgaWYgKG9yZGVyID09PSAnWlhaJykge1xuICAgIC8vIGF4aXNBbmdsZShbMCwgMCwgMV0sIM+GKSAqIGF4aXNBbmdsZShbMSwgMCwgMF0sIM64KSAqIGF4aXNBbmdsZShbMCwgMCwgMV0sIM+IKVxuICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKFxuICAgICAgY1ggKiBjWSAqIGNaIC0gc1ggKiBzWiAqIGNZLFxuICAgICAgc1ggKiBzWSAqIHNaICsgc1kgKiBjWCAqIGNaLFxuICAgICAgc1ggKiBzWSAqIGNaIC0gc1kgKiBzWiAqIGNYLFxuICAgICAgc1ggKiBjWSAqIGNaICsgc1ogKiBjWCAqIGNZKTtcbiAgfVxuXG4gIGlmIChvcmRlciA9PT0gJ1lYWScpIHtcbiAgICAvLyBheGlzQW5nbGUoWzAsIDEsIDBdLCDPhikgKiBheGlzQW5nbGUoWzEsIDAsIDBdLCDOuCkgKiBheGlzQW5nbGUoWzAsIDEsIDBdLCDPiClcbiAgICByZXR1cm4gbmV3UXVhdGVybmlvbihcbiAgICAgIGNYICogY1kgKiBjWiAtIHNYICogc1ogKiBjWSxcbiAgICAgIHNYICogc1kgKiBzWiArIHNZICogY1ggKiBjWixcbiAgICAgIHNYICogY1kgKiBjWiArIHNaICogY1ggKiBjWSxcbiAgICAgIHNZICogc1ogKiBjWCAtIHNYICogc1kgKiBjWik7XG4gIH1cblxuICBpZiAob3JkZXIgPT09ICdZWlknKSB7XG4gICAgLy8gYXhpc0FuZ2xlKFswLCAxLCAwXSwgz4YpICogYXhpc0FuZ2xlKFswLCAwLCAxXSwgzrgpICogYXhpc0FuZ2xlKFswLCAxLCAwXSwgz4gpXG4gICAgcmV0dXJuIG5ld1F1YXRlcm5pb24oXG4gICAgICBjWCAqIGNZICogY1ogLSBzWCAqIHNaICogY1ksXG4gICAgICBzWCAqIHNZICogY1ogLSBzWSAqIHNaICogY1gsXG4gICAgICBzWCAqIGNZICogY1ogKyBzWiAqIGNYICogY1ksXG4gICAgICBzWCAqIHNZICogc1ogKyBzWSAqIGNYICogY1opO1xuICB9XG5cbiAgaWYgKG9yZGVyID09PSAnWFlYJykge1xuICAgIC8vIGF4aXNBbmdsZShbMSwgMCwgMF0sIM+GKSAqIGF4aXNBbmdsZShbMCwgMSwgMF0sIM64KSAqIGF4aXNBbmdsZShbMSwgMCwgMF0sIM+IKVxuICAgIHJldHVybiBuZXdRdWF0ZXJuaW9uKFxuICAgICAgY1ggKiBjWSAqIGNaIC0gc1ggKiBzWiAqIGNZLFxuICAgICAgc1ggKiBjWSAqIGNaICsgc1ogKiBjWCAqIGNZLFxuICAgICAgc1ggKiBzWSAqIHNaICsgc1kgKiBjWCAqIGNaLFxuICAgICAgc1ggKiBzWSAqIGNaIC0gc1kgKiBzWiAqIGNYKTtcbiAgfVxuXG4gIGlmIChvcmRlciA9PT0gJ1haWCcpIHtcbiAgICAvLyBheGlzQW5nbGUoWzEsIDAsIDBdLCDPhikgKiBheGlzQW5nbGUoWzAsIDAsIDFdLCDOuCkgKiBheGlzQW5nbGUoWzEsIDAsIDBdLCDPiClcbiAgICByZXR1cm4gbmV3UXVhdGVybmlvbihcbiAgICAgIGNYICogY1kgKiBjWiAtIHNYICogc1ogKiBjWSxcbiAgICAgIHNYICogY1kgKiBjWiArIHNaICogY1ggKiBjWSxcbiAgICAgIHNZICogc1ogKiBjWCAtIHNYICogc1kgKiBjWixcbiAgICAgIHNYICogc1kgKiBzWiArIHNZICogY1ggKiBjWik7XG4gIH1cbiAgcmV0dXJuIG51bGw7XG59O1xuXG4vKipcbiAqIENyZWF0ZXMgYSBxdWF0ZXJuaW9uIGJ5IGEgcm90YXRpb24gbWF0cml4XG4gKlxuICogQHBhcmFtIHtBcnJheX0gbWF0cml4XG4gKiBAcmV0dXJucyB7UXVhdGVybmlvbn1cbiAqL1xuUXVhdGVybmlvblsnZnJvbU1hdHJpeCddID0gZnVuY3Rpb24gKG1hdHJpeCkge1xuXG4gIGxldCBtMDAsIG0wMSwgbTAyLCBtMTAsIG0xMSwgbTEyLCBtMjAsIG0yMSwgbTIyO1xuXG4gIGlmIChtYXRyaXgubGVuZ3RoID09PSA5KSB7XG4gICAgbTAwID0gbWF0cml4WzBdO1xuICAgIG0wMSA9IG1hdHJpeFsxXTtcbiAgICBtMDIgPSBtYXRyaXhbMl07XG5cbiAgICBtMTAgPSBtYXRyaXhbM107XG4gICAgbTExID0gbWF0cml4WzRdO1xuICAgIG0xMiA9IG1hdHJpeFs1XTtcblxuICAgIG0yMCA9IG1hdHJpeFs2XTtcbiAgICBtMjEgPSBtYXRyaXhbN107XG4gICAgbTIyID0gbWF0cml4WzhdO1xuXG4gIH0gZWxzZSB7XG4gICAgbTAwID0gbWF0cml4WzBdWzBdO1xuICAgIG0wMSA9IG1hdHJpeFswXVsxXTtcbiAgICBtMDIgPSBtYXRyaXhbMF1bMl07XG5cbiAgICBtMTAgPSBtYXRyaXhbMV1bMF07XG4gICAgbTExID0gbWF0cml4WzFdWzFdO1xuICAgIG0xMiA9IG1hdHJpeFsxXVsyXTtcblxuICAgIG0yMCA9IG1hdHJpeFsyXVswXTtcbiAgICBtMjEgPSBtYXRyaXhbMl1bMV07XG4gICAgbTIyID0gbWF0cml4WzJdWzJdO1xuICB9XG5cbiAgY29uc3QgdHIgPSBtMDAgKyBtMTEgKyBtMjI7IC8vIDIgKiB3ID0gc3FydCgxICsgdHIpXG5cbiAgLy8gQ2hvb3NlIHRoZSBlbGVtZW50IHdpdGggdGhlIGJpZ2dlc3QgdmFsdWUgb24gdGhlIGRpYWdvbmFsXG5cbiAgaWYgKHRyID4gMCkgeyAvLyBpZiB0cmFjZSBpcyBwb3NpdGl2ZSB0aGVuIFwid1wiIGlzIGJpZ2dlc3QgY29tcG9uZW50XG4gICAgLy8gfFF8ID0gMiAqIHNxcnQoMSArIHRyKSA9IDR3XG4gICAgcmV0dXJuIG5ld05vcm1hbGl6ZWQoXG4gICAgICB0ciArIDEuMCxcbiAgICAgIG0yMSAtIG0xMixcbiAgICAgIG0wMiAtIG0yMCxcbiAgICAgIG0xMCAtIG0wMSk7XG4gIH0gZWxzZSBpZiAobTAwID4gbTExICYmIG0wMCA+IG0yMikge1xuICAgIC8vIHxRfCA9IDIgKiBzcXJ0KDEuMCArIG0wMCAtIG0xMSAtIG0yMikgPSA0eFxuICAgIHJldHVybiBuZXdOb3JtYWxpemVkKFxuICAgICAgbTIxIC0gbTEyLFxuICAgICAgMS4wICsgbTAwIC0gbTExIC0gbTIyLFxuICAgICAgbTAxICsgbTEwLFxuICAgICAgbTAyICsgbTIwKTtcbiAgfSBlbHNlIGlmIChtMTEgPiBtMjIpIHtcbiAgICAvLyB8UXwgPSAyICogc3FydCgxLjAgKyBtMTEgLSBtMDAgLSBtMjIpID0gNHlcbiAgICByZXR1cm4gbmV3Tm9ybWFsaXplZChcbiAgICAgIG0wMiAtIG0yMCxcbiAgICAgIG0wMSArIG0xMCxcbiAgICAgIDEuMCArIG0xMSAtIG0wMCAtIG0yMixcbiAgICAgIG0xMiArIG0yMSk7XG4gIH0gZWxzZSB7XG4gICAgLy8gfFF8ID0gMiAqIHNxcnQoMS4wICsgbTIyIC0gbTAwIC0gbTExKSA9IDR6XG4gICAgcmV0dXJuIG5ld05vcm1hbGl6ZWQoXG4gICAgICBtMTAgLSBtMDEsXG4gICAgICBtMDIgKyBtMjAsXG4gICAgICBtMTIgKyBtMjEsXG4gICAgICAxLjAgKyBtMjIgLSBtMDAgLSBtMTEpO1xuICB9XG59O1xuZXhwb3J0IHtcbiAgUXVhdGVybmlvbiBhcyBkZWZhdWx0LCBRdWF0ZXJuaW9uXG59O1xuIiwiZXhwb3J0IGNvbnN0IGFsaWFzZXMgPSB7XG4gICAgXCJtUGVsdmlzXCI6IFwiaGlwXCIsXG4gICAgXCJtU3BpbmUxXCI6IFwibVNwaW5lMVwiLFxuICAgIFwibVNwaW5lMlwiOiBcIm1TcGluZTJcIixcbiAgICBcIm1Ub3Jzb1wiOiBcImFiZG9tZW5cIixcbiAgICBcIm1TcGluZTNcIjogXCJtU3BpbmUzXCIsXG4gICAgXCJtU3BpbmU0XCI6IFwibVNwaW5lNFwiLFxuICAgIFwibUNoZXN0XCI6IFwiY2hlc3RcIixcbiAgICBcIm1OZWNrXCI6IFwibmVja1wiLFxuICAgIFwibUhlYWRcIjogXCJoZWFkXCIsXG4gICAgXCJtU2t1bGxcIjogXCJmaWd1cmVIYWlyXCIsXG4gICAgXCJtRXllUmlnaHRcIjogXCJtRXllUmlnaHRcIixcbiAgICBcIm1FeWVMZWZ0XCI6IFwibUV5ZUxlZnRcIixcbiAgICBcIm1GYWNlUm9vdFwiOiBcIm1GYWNlUm9vdFwiLFxuICAgIFwibUZhY2VFeWVBbHRSaWdodFwiOiBcIm1GYWNlRXllQWx0UmlnaHRcIixcbiAgICBcIm1GYWNlRXllQWx0TGVmdFwiOiBcIm1GYWNlRXllQWx0TGVmdFwiLFxuICAgIFwibUZhY2VGb3JlaGVhZExlZnRcIjogXCJtRmFjZUZvcmVoZWFkTGVmdFwiLFxuICAgIFwibUZhY2VGb3JlaGVhZFJpZ2h0XCI6IFwibUZhY2VGb3JlaGVhZFJpZ2h0XCIsXG4gICAgXCJtRmFjZUV5ZWJyb3dPdXRlckxlZnRcIjogXCJtRmFjZUV5ZWJyb3dPdXRlckxlZnRcIixcbiAgICBcIm1GYWNlRXllYnJvd0NlbnRlckxlZnRcIjogXCJtRmFjZUV5ZWJyb3dDZW50ZXJMZWZ0XCIsXG4gICAgXCJtRmFjZUV5ZWJyb3dJbm5lckxlZnRcIjogXCJtRmFjZUV5ZWJyb3dJbm5lckxlZnRcIixcbiAgICBcIm1GYWNlRXllYnJvd091dGVyUmlnaHRcIjogXCJtRmFjZUV5ZWJyb3dPdXRlclJpZ2h0XCIsXG4gICAgXCJtRmFjZUV5ZWJyb3dDZW50ZXJSaWdodFwiOiBcIm1GYWNlRXllYnJvd0NlbnRlclJpZ2h0XCIsXG4gICAgXCJtRmFjZUV5ZWJyb3dJbm5lclJpZ2h0XCI6IFwibUZhY2VFeWVicm93SW5uZXJSaWdodFwiLFxuICAgIFwibUZhY2VFeWVMaWRVcHBlckxlZnRcIjogXCJtRmFjZUV5ZUxpZFVwcGVyTGVmdFwiLFxuICAgIFwibUZhY2VFeWVMaWRMb3dlckxlZnRcIjogXCJtRmFjZUV5ZUxpZExvd2VyTGVmdFwiLFxuICAgIFwibUZhY2VFeWVMaWRVcHBlclJpZ2h0XCI6IFwibUZhY2VFeWVMaWRVcHBlclJpZ2h0XCIsXG4gICAgXCJtRmFjZUV5ZUxpZExvd2VyUmlnaHRcIjogXCJtRmFjZUV5ZUxpZExvd2VyUmlnaHRcIixcbiAgICBcIm1GYWNlRWFyMUxlZnRcIjogXCJtRmFjZUVhcjFMZWZ0XCIsXG4gICAgXCJtRmFjZUVhcjJMZWZ0XCI6IFwibUZhY2VFYXIyTGVmdFwiLFxuICAgIFwibUZhY2VFYXIxUmlnaHRcIjogXCJtRmFjZUVhcjFSaWdodFwiLFxuICAgIFwibUZhY2VFYXIyUmlnaHRcIjogXCJtRmFjZUVhcjJSaWdodFwiLFxuICAgIFwibUZhY2VOb3NlTGVmdFwiOiBcIm1GYWNlTm9zZUxlZnRcIixcbiAgICBcIm1GYWNlTm9zZUNlbnRlclwiOiBcIm1GYWNlTm9zZUNlbnRlclwiLFxuICAgIFwibUZhY2VOb3NlUmlnaHRcIjogXCJtRmFjZU5vc2VSaWdodFwiLFxuICAgIFwibUZhY2VDaGVla0xvd2VyTGVmdFwiOiBcIm1GYWNlQ2hlZWtMb3dlckxlZnRcIixcbiAgICBcIm1GYWNlQ2hlZWtVcHBlckxlZnRcIjogXCJtRmFjZUNoZWVrVXBwZXJMZWZ0XCIsXG4gICAgXCJtRmFjZUNoZWVrTG93ZXJSaWdodFwiOiBcIm1GYWNlQ2hlZWtMb3dlclJpZ2h0XCIsXG4gICAgXCJtRmFjZUNoZWVrVXBwZXJSaWdodFwiOiBcIm1GYWNlQ2hlZWtVcHBlclJpZ2h0XCIsXG4gICAgXCJtRmFjZUphd1wiOiBcIm1GYWNlSmF3XCIsXG4gICAgXCJtRmFjZUNoaW5cIjogXCJtRmFjZUNoaW5cIixcbiAgICBcIm1GYWNlVGVldGhMb3dlclwiOiBcIm1GYWNlVGVldGhMb3dlclwiLFxuICAgIFwibUZhY2VMaXBMb3dlckxlZnRcIjogXCJtRmFjZUxpcExvd2VyTGVmdFwiLFxuICAgIFwibUZhY2VMaXBMb3dlclJpZ2h0XCI6IFwibUZhY2VMaXBMb3dlclJpZ2h0XCIsXG4gICAgXCJtRmFjZUxpcExvd2VyQ2VudGVyXCI6IFwibUZhY2VMaXBMb3dlckNlbnRlclwiLFxuICAgIFwibUZhY2VUb25ndWVCYXNlXCI6IFwibUZhY2VUb25ndWVCYXNlXCIsXG4gICAgXCJtRmFjZVRvbmd1ZVRpcFwiOiBcIm1GYWNlVG9uZ3VlVGlwXCIsXG4gICAgXCJtRmFjZUphd1NoYXBlclwiOiBcIm1GYWNlSmF3U2hhcGVyXCIsXG4gICAgXCJtRmFjZUZvcmVoZWFkQ2VudGVyXCI6IFwibUZhY2VGb3JlaGVhZENlbnRlclwiLFxuICAgIFwibUZhY2VOb3NlQmFzZVwiOiBcIm1GYWNlTm9zZUJhc2VcIixcbiAgICBcIm1GYWNlVGVldGhVcHBlclwiOiBcIm1GYWNlVGVldGhVcHBlclwiLFxuICAgIFwibUZhY2VMaXBVcHBlckxlZnRcIjogXCJtRmFjZUxpcFVwcGVyTGVmdFwiLFxuICAgIFwibUZhY2VMaXBVcHBlclJpZ2h0XCI6IFwibUZhY2VMaXBVcHBlclJpZ2h0XCIsXG4gICAgXCJtRmFjZUxpcENvcm5lckxlZnRcIjogXCJtRmFjZUxpcENvcm5lckxlZnRcIixcbiAgICBcIm1GYWNlTGlwQ29ybmVyUmlnaHRcIjogXCJtRmFjZUxpcENvcm5lclJpZ2h0XCIsXG4gICAgXCJtRmFjZUxpcFVwcGVyQ2VudGVyXCI6IFwibUZhY2VMaXBVcHBlckNlbnRlclwiLFxuICAgIFwibUZhY2VFeWVjb3JuZXJJbm5lckxlZnRcIjogXCJtRmFjZUV5ZWNvcm5lcklubmVyTGVmdFwiLFxuICAgIFwibUZhY2VFeWVjb3JuZXJJbm5lclJpZ2h0XCI6IFwibUZhY2VFeWVjb3JuZXJJbm5lclJpZ2h0XCIsXG4gICAgXCJtRmFjZU5vc2VCcmlkZ2VcIjogXCJtRmFjZU5vc2VCcmlkZ2VcIixcbiAgICBcIm1Db2xsYXJMZWZ0XCI6IFwibENvbGxhclwiLFxuICAgIFwibVNob3VsZGVyTGVmdFwiOiBcImxTaGxkclwiLFxuICAgIFwibUVsYm93TGVmdFwiOiBcImxGb3JlQXJtXCIsXG4gICAgXCJtV3Jpc3RMZWZ0XCI6IFwibEhhbmRcIixcbiAgICBcIm1IYW5kTWlkZGxlMUxlZnRcIjogXCJtSGFuZE1pZGRsZTFMZWZ0XCIsXG4gICAgXCJtSGFuZE1pZGRsZTJMZWZ0XCI6IFwibUhhbmRNaWRkbGUyTGVmdFwiLFxuICAgIFwibUhhbmRNaWRkbGUzTGVmdFwiOiBcIm1IYW5kTWlkZGxlM0xlZnRcIixcbiAgICBcIm1IYW5kSW5kZXgxTGVmdFwiOiBcIm1IYW5kSW5kZXgxTGVmdFwiLFxuICAgIFwibUhhbmRJbmRleDJMZWZ0XCI6IFwibUhhbmRJbmRleDJMZWZ0XCIsXG4gICAgXCJtSGFuZEluZGV4M0xlZnRcIjogXCJtSGFuZEluZGV4M0xlZnRcIixcbiAgICBcIm1IYW5kUmluZzFMZWZ0XCI6IFwibUhhbmRSaW5nMUxlZnRcIixcbiAgICBcIm1IYW5kUmluZzJMZWZ0XCI6IFwibUhhbmRSaW5nMkxlZnRcIixcbiAgICBcIm1IYW5kUmluZzNMZWZ0XCI6IFwibUhhbmRSaW5nM0xlZnRcIixcbiAgICBcIm1IYW5kUGlua3kxTGVmdFwiOiBcIm1IYW5kUGlua3kxTGVmdFwiLFxuICAgIFwibUhhbmRQaW5reTJMZWZ0XCI6IFwibUhhbmRQaW5reTJMZWZ0XCIsXG4gICAgXCJtSGFuZFBpbmt5M0xlZnRcIjogXCJtSGFuZFBpbmt5M0xlZnRcIixcbiAgICBcIm1IYW5kVGh1bWIxTGVmdFwiOiBcIm1IYW5kVGh1bWIxTGVmdFwiLFxuICAgIFwibUhhbmRUaHVtYjJMZWZ0XCI6IFwibUhhbmRUaHVtYjJMZWZ0XCIsXG4gICAgXCJtSGFuZFRodW1iM0xlZnRcIjogXCJtSGFuZFRodW1iM0xlZnRcIixcbiAgICBcIm1Db2xsYXJSaWdodFwiOiBcInJDb2xsYXJcIixcbiAgICBcIm1TaG91bGRlclJpZ2h0XCI6IFwiclNobGRyXCIsXG4gICAgXCJtRWxib3dSaWdodFwiOiBcInJGb3JlQXJtXCIsXG4gICAgXCJtV3Jpc3RSaWdodFwiOiBcInJIYW5kXCIsXG4gICAgXCJtSGFuZE1pZGRsZTFSaWdodFwiOiBcIm1IYW5kTWlkZGxlMVJpZ2h0XCIsXG4gICAgXCJtSGFuZE1pZGRsZTJSaWdodFwiOiBcIm1IYW5kTWlkZGxlMlJpZ2h0XCIsXG4gICAgXCJtSGFuZE1pZGRsZTNSaWdodFwiOiBcIm1IYW5kTWlkZGxlM1JpZ2h0XCIsXG4gICAgXCJtSGFuZEluZGV4MVJpZ2h0XCI6IFwibUhhbmRJbmRleDFSaWdodFwiLFxuICAgIFwibUhhbmRJbmRleDJSaWdodFwiOiBcIm1IYW5kSW5kZXgyUmlnaHRcIixcbiAgICBcIm1IYW5kSW5kZXgzUmlnaHRcIjogXCJtSGFuZEluZGV4M1JpZ2h0XCIsXG4gICAgXCJtSGFuZFJpbmcxUmlnaHRcIjogXCJtSGFuZFJpbmcxUmlnaHRcIixcbiAgICBcIm1IYW5kUmluZzJSaWdodFwiOiBcIm1IYW5kUmluZzJSaWdodFwiLFxuICAgIFwibUhhbmRSaW5nM1JpZ2h0XCI6IFwibUhhbmRSaW5nM1JpZ2h0XCIsXG4gICAgXCJtSGFuZFBpbmt5MVJpZ2h0XCI6IFwibUhhbmRQaW5reTFSaWdodFwiLFxuICAgIFwibUhhbmRQaW5reTJSaWdodFwiOiBcIm1IYW5kUGlua3kyUmlnaHRcIixcbiAgICBcIm1IYW5kUGlua3kzUmlnaHRcIjogXCJtSGFuZFBpbmt5M1JpZ2h0XCIsXG4gICAgXCJtSGFuZFRodW1iMVJpZ2h0XCI6IFwibUhhbmRUaHVtYjFSaWdodFwiLFxuICAgIFwibUhhbmRUaHVtYjJSaWdodFwiOiBcIm1IYW5kVGh1bWIyUmlnaHRcIixcbiAgICBcIm1IYW5kVGh1bWIzUmlnaHRcIjogXCJtSGFuZFRodW1iM1JpZ2h0XCIsXG4gICAgXCJtV2luZ3NSb290XCI6IFwibVdpbmdzUm9vdFwiLFxuICAgIFwibVdpbmcxTGVmdFwiOiBcIm1XaW5nMUxlZnRcIixcbiAgICBcIm1XaW5nMkxlZnRcIjogXCJtV2luZzJMZWZ0XCIsXG4gICAgXCJtV2luZzNMZWZ0XCI6IFwibVdpbmczTGVmdFwiLFxuICAgIFwibVdpbmc0TGVmdFwiOiBcIm1XaW5nNExlZnRcIixcbiAgICBcIm1XaW5nNEZhbkxlZnRcIjogXCJtV2luZzRGYW5MZWZ0XCIsXG4gICAgXCJtV2luZzFSaWdodFwiOiBcIm1XaW5nMVJpZ2h0XCIsXG4gICAgXCJtV2luZzJSaWdodFwiOiBcIm1XaW5nMlJpZ2h0XCIsXG4gICAgXCJtV2luZzNSaWdodFwiOiBcIm1XaW5nM1JpZ2h0XCIsXG4gICAgXCJtV2luZzRSaWdodFwiOiBcIm1XaW5nNFJpZ2h0XCIsXG4gICAgXCJtV2luZzRGYW5SaWdodFwiOiBcIm1XaW5nNEZhblJpZ2h0XCIsXG4gICAgXCJtSGlwUmlnaHRcIjogXCJyVGhpZ2hcIixcbiAgICBcIm1LbmVlUmlnaHRcIjogXCJyU2hpblwiLFxuICAgIFwibUFua2xlUmlnaHRcIjogXCJyRm9vdFwiLFxuICAgIFwibUZvb3RSaWdodFwiOiBcIm1Gb290UmlnaHRcIixcbiAgICBcIm1Ub2VSaWdodFwiOiBcIm1Ub2VSaWdodFwiLFxuICAgIFwibUhpcExlZnRcIjogXCJsVGhpZ2hcIixcbiAgICBcIm1LbmVlTGVmdFwiOiBcImxTaGluXCIsXG4gICAgXCJtQW5rbGVMZWZ0XCI6IFwibEZvb3RcIixcbiAgICBcIm1Gb290TGVmdFwiOiBcIm1Gb290TGVmdFwiLFxuICAgIFwibVRvZUxlZnRcIjogXCJtVG9lTGVmdFwiLFxuICAgIFwibVRhaWwxXCI6IFwibVRhaWwxXCIsXG4gICAgXCJtVGFpbDJcIjogXCJtVGFpbDJcIixcbiAgICBcIm1UYWlsM1wiOiBcIm1UYWlsM1wiLFxuICAgIFwibVRhaWw0XCI6IFwibVRhaWw0XCIsXG4gICAgXCJtVGFpbDVcIjogXCJtVGFpbDVcIixcbiAgICBcIm1UYWlsNlwiOiBcIm1UYWlsNlwiLFxuICAgIFwibUdyb2luXCI6IFwibUdyb2luXCIsXG4gICAgXCJtSGluZExpbWJzUm9vdFwiOiBcIm1IaW5kTGltYnNSb290XCIsXG4gICAgXCJtSGluZExpbWIxTGVmdFwiOiBcIm1IaW5kTGltYjFMZWZ0XCIsXG4gICAgXCJtSGluZExpbWIyTGVmdFwiOiBcIm1IaW5kTGltYjJMZWZ0XCIsXG4gICAgXCJtSGluZExpbWIzTGVmdFwiOiBcIm1IaW5kTGltYjNMZWZ0XCIsXG4gICAgXCJtSGluZExpbWI0TGVmdFwiOiBcIm1IaW5kTGltYjRMZWZ0XCIsXG4gICAgXCJtSGluZExpbWIxUmlnaHRcIjogXCJtSGluZExpbWIxUmlnaHRcIixcbiAgICBcIm1IaW5kTGltYjJSaWdodFwiOiBcIm1IaW5kTGltYjJSaWdodFwiLFxuICAgIFwibUhpbmRMaW1iM1JpZ2h0XCI6IFwibUhpbmRMaW1iM1JpZ2h0XCIsXG4gICAgXCJtSGluZExpbWI0UmlnaHRcIjogXCJtSGluZExpbWI0UmlnaHRcIlxufTtcbiIsImltcG9ydCB7IHRvUXVhdGVybmlvbiwgbGVycFZhbHVlcywgZ2V0VW5pZm9ybVRpbWVzLCBjbGlwVGltZXNUb0Nsb3Nlc3RCVkhUaW1lLCBsZXJwVmVjdG9yLCBsZXJwUXVhdGVybmlvbiwgcXVhdGVybmlvblRvRXVsZXJzLCBmbG9hdFRvU3RyaW5nIH0gZnJvbSBcIi4vdXRpbHNcIjtcbmltcG9ydCB7IGhpZXJhcmNoeSB9IGZyb20gXCIuL2hpZXJhcmNoeVwiO1xuaW1wb3J0IHsgYWxpYXNlcyB9IGZyb20gXCIuL2FsaWFzZXNcIjtcbmZ1bmN0aW9uIG9mZnNldFRvU3RyaW5nKG9mZnNldCwgZGlnaXRzKSB7XG4gICAgcmV0dXJuIFwiT0ZGU0VUIFwiICsgZmxvYXRUb1N0cmluZyhvZmZzZXQueCwgZGlnaXRzKSArIFwiIFwiICsgZmxvYXRUb1N0cmluZyhvZmZzZXQueSwgZGlnaXRzKSArIFwiIFwiICsgZmxvYXRUb1N0cmluZyhvZmZzZXQueiwgZGlnaXRzKTtcbn1cbmZ1bmN0aW9uIGFwcGVuZE5vZGUoam9pbnQsIHRhYnMpIHtcbiAgICBsZXQgcmVzdWx0ID0gXCJcIjtcbiAgICBjb25zdCBib25lVHlwZSA9IChqb2ludC5idmhOYW1lID09PSBcImhpcFwiKSA/IFwiUk9PVFwiIDogXCJKT0lOVFwiO1xuICAgIGNvbnN0IGNoYW5uZWxzID0gKGpvaW50LmJ2aE5hbWUgPT09IFwiaGlwXCIpID8gXCJDSEFOTkVMUyA2IFhwb3NpdGlvbiBZcG9zaXRpb24gWnBvc2l0aW9uIFhyb3RhdGlvbiBZcm90YXRpb24gWnJvdGF0aW9uXCIgOiBcIkNIQU5ORUxTIDMgWHJvdGF0aW9uIFlyb3RhdGlvbiBacm90YXRpb25cIjtcbiAgICBjb25zdCBvZmZzZXQgPSAoam9pbnQuYnZoTmFtZSA9PT0gXCJoaXBcIikgPyBvZmZzZXRUb1N0cmluZyhqb2ludC5vZmZzZXQsIDYpIDogb2Zmc2V0VG9TdHJpbmcoam9pbnQub2Zmc2V0LCA0KTtcbiAgICBpZiAoam9pbnQuYnZoTmFtZSAhPSBcImVuZFwiKSB7XG4gICAgICAgIHJlc3VsdCArPSB0YWJzICsgYm9uZVR5cGUgKyBcIiBcIiArIGpvaW50LmJ2aE5hbWUgKyBcIlxcblwiICsgdGFicyArIFwie1xcblwiO1xuICAgIH1cbiAgICBlbHNlIHtcbiAgICAgICAgcmVzdWx0ICs9IHRhYnMgKyBcIkVuZCBTaXRlXCIgKyBcIlxcblwiICsgdGFicyArIFwie1xcblwiO1xuICAgIH1cbiAgICByZXN1bHQgKz0gdGFicyArIFwiXFx0XCIgKyBvZmZzZXQgKyBcIlxcblwiO1xuICAgIGlmIChqb2ludC5idmhOYW1lICE9IFwiZW5kXCIpIHtcbiAgICAgICAgcmVzdWx0ICs9IHRhYnMgKyBcIlxcdFwiICsgY2hhbm5lbHMgKyBcIlxcblwiO1xuICAgIH1cbiAgICBpZiAoam9pbnQuY2hpbGRyZW4pIHtcbiAgICAgICAgam9pbnQuY2hpbGRyZW4uZm9yRWFjaCgoaXRlbSkgPT4geyByZXN1bHQgKz0gYXBwZW5kTm9kZShpdGVtLCB0YWJzICsgXCJcXHRcIik7IH0pO1xuICAgIH1cbiAgICByZXN1bHQgKz0gdGFicyArIFwifVxcblwiO1xuICAgIHJldHVybiByZXN1bHQ7XG59XG5mdW5jdGlvbiBjb250YWluc05hbWVzKG5vZGUsIGJ2aE5hbWVzKSB7XG4gICAgaWYgKGJ2aE5hbWVzLmluY2x1ZGVzKG5vZGUuYnZoTmFtZSkpIHtcbiAgICAgICAgcmV0dXJuIHRydWU7XG4gICAgfVxuICAgIGlmICghbm9kZS5jaGlsZHJlbikge1xuICAgICAgICByZXR1cm4gZmFsc2U7XG4gICAgfVxuICAgIHJldHVybiAhIW5vZGUuY2hpbGRyZW4ubWFwKChpdGVtKSA9PiBjb250YWluc05hbWVzKGl0ZW0sIGJ2aE5hbWVzKSkuZmluZCgoaXRlbSkgPT4gISFpdGVtKTtcbn1cbmZ1bmN0aW9uIGNvbGxlY3ROb2Rlcyhub2RlLCBidmhOYW1lcykge1xuICAgIGNvbnN0IHJlc3VsdCA9IHt9O1xuICAgIGlmIChjb250YWluc05hbWVzKG5vZGUsIGJ2aE5hbWVzKSkge1xuICAgICAgICByZXN1bHQuYnZoTmFtZSA9IG5vZGUuYnZoTmFtZTtcbiAgICB9XG4gICAgZWxzZSB7XG4gICAgICAgIHJlc3VsdC5leGNsdWRlID0gdHJ1ZTtcbiAgICAgICAgcmV0dXJuIHJlc3VsdDtcbiAgICB9XG4gICAgaWYgKG5vZGUuY2hpbGRyZW4gJiYgISFub2RlLmNoaWxkcmVuLm1hcCgoaXRlbSkgPT4gY29udGFpbnNOYW1lcyhpdGVtLCBidmhOYW1lcykpLmZpbmQoKGl0ZW0pID0+ICEhaXRlbSkpIHtcbiAgICAgICAgcmVzdWx0LmNoaWxkcmVuID0gbm9kZS5jaGlsZHJlbi5tYXAoKGl0ZW0pID0+IGNvbGxlY3ROb2RlcyhpdGVtLCBidmhOYW1lcykpLmZpbHRlcigoaXRlbSkgPT4gIWl0ZW0uZXhjbHVkZSk7XG4gICAgfVxuICAgIGVsc2Uge1xuICAgICAgICByZXN1bHQuY2hpbGRyZW4gPSBbXTtcbiAgICB9XG4gICAgaWYgKHJlc3VsdC5jaGlsZHJlbi5sZW5ndGggPiAwKSB7XG4gICAgICAgIHJldHVybiByZXN1bHQ7XG4gICAgfVxuICAgIHJlc3VsdC5jaGlsZHJlbi5wdXNoKHsgYnZoTmFtZTogXCJlbmRcIiB9KTtcbiAgICByZXR1cm4gcmVzdWx0O1xufVxuZnVuY3Rpb24gc3ViVHJlZShqb2ludHMpIHtcbiAgICBjb25zdCBuYW1lcyA9IGpvaW50cy5tYXAoaXRlbSA9PiBpdGVtLmpvaW50X25hbWUpO1xuICAgIGNvbnN0IGJ2aE5hbWVzID0gbmFtZXMubWFwKGl0ZW0gPT4gYWxpYXNlc1tpdGVtXSB8fCBpdGVtKTtcbiAgICByZXR1cm4gY29sbGVjdE5vZGVzKGhpZXJhcmNoeSwgYnZoTmFtZXMpO1xufVxuZXhwb3J0IGZ1bmN0aW9uIHZpc2l0Tm9kZShub2RlLCB2aXNpdG9yLCBjaGlsZHJlbkZpcnN0ID0gZmFsc2UpIHtcbiAgICBpZiAobm9kZS5jaGlsZHJlbiAmJiBjaGlsZHJlbkZpcnN0KSB7XG4gICAgICAgIG5vZGUuY2hpbGRyZW4udG9SZXZlcnNlZCgpLmZvckVhY2goKGl0ZW0pID0+IHZpc2l0Tm9kZShpdGVtLCB2aXNpdG9yLCB0cnVlKSk7XG4gICAgfVxuICAgIHZpc2l0b3Iobm9kZSk7XG4gICAgaWYgKG5vZGUuY2hpbGRyZW4gJiYgIWNoaWxkcmVuRmlyc3QpIHtcbiAgICAgICAgbm9kZS5jaGlsZHJlbi50b1JldmVyc2VkKCkuZm9yRWFjaCgoaXRlbSkgPT4gdmlzaXROb2RlKGl0ZW0sIHZpc2l0b3IsIGZhbHNlKSk7XG4gICAgfVxufVxuZnVuY3Rpb24gZXh0cmFjdEZyYW1lc0xlbmd0aChhbmltSm9pbnRzKSB7XG4gICAgdmFyIF9hLCBfYjtcbiAgICBjb25zdCBqb2ludCA9IGFuaW1Kb2ludHMuZmluZCgoaXRlbSkgPT4geyB2YXIgX2EsIF9iOyByZXR1cm4gKChfYSA9IGl0ZW0ucG9zaXRpb25fa2V5cykgPT09IG51bGwgfHwgX2EgPT09IHZvaWQgMCA/IHZvaWQgMCA6IF9hLmxlbmd0aCkgfHwgKChfYiA9IGl0ZW0ucm90YXRpb25fa2V5cykgPT09IG51bGwgfHwgX2IgPT09IHZvaWQgMCA/IHZvaWQgMCA6IF9iLmxlbmd0aCk7IH0pO1xuICAgIHJldHVybiAoKF9hID0gam9pbnQgPT09IG51bGwgfHwgam9pbnQgPT09IHZvaWQgMCA/IHZvaWQgMCA6IGpvaW50LnBvc2l0aW9uX2tleXMpID09PSBudWxsIHx8IF9hID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfYS5sZW5ndGgpIHx8ICgoX2IgPSBqb2ludCA9PT0gbnVsbCB8fCBqb2ludCA9PT0gdm9pZCAwID8gdm9pZCAwIDogam9pbnQucm90YXRpb25fa2V5cykgPT09IG51bGwgfHwgX2IgPT09IHZvaWQgMCA/IHZvaWQgMCA6IF9iLmxlbmd0aCk7XG59XG5mdW5jdGlvbiBleHRyYWN0VGltZXMoYW5pbUpvaW50cykge1xuICAgIHZhciBfYTtcbiAgICBjb25zdCBqb2ludCA9IGFuaW1Kb2ludHMuZmluZCgoaXRlbSkgPT4geyB2YXIgX2EsIF9iOyByZXR1cm4gKChfYSA9IGl0ZW0ucG9zaXRpb25fa2V5cykgPT09IG51bGwgfHwgX2EgPT09IHZvaWQgMCA/IHZvaWQgMCA6IF9hLmxlbmd0aCkgfHwgKChfYiA9IGl0ZW0ucm90YXRpb25fa2V5cykgPT09IG51bGwgfHwgX2IgPT09IHZvaWQgMCA/IHZvaWQgMCA6IF9iLmxlbmd0aCk7IH0pO1xuICAgIGNvbnN0IHRpbWVIb2xkZXJzID0gKChfYSA9IGpvaW50ID09PSBudWxsIHx8IGpvaW50ID09PSB2b2lkIDAgPyB2b2lkIDAgOiBqb2ludC5wb3NpdGlvbl9rZXlzKSA9PT0gbnVsbCB8fCBfYSA9PT0gdm9pZCAwID8gdm9pZCAwIDogX2EubGVuZ3RoKSA/IGpvaW50LnBvc2l0aW9uX2tleXMgOiBqb2ludCA9PT0gbnVsbCB8fCBqb2ludCA9PT0gdm9pZCAwID8gdm9pZCAwIDogam9pbnQucm90YXRpb25fa2V5cztcbiAgICByZXR1cm4gKHRpbWVIb2xkZXJzIHx8IFtdKS5tYXAoKGl0ZW0pID0+IGl0ZW0udGltZSk7XG59XG5mdW5jdGlvbiBmaWxsQ2hhbm5lbHMobm9kZSwgam9pbnQpIHtcbiAgICB2YXIgX2E7XG4gICAgaWYgKG5vZGUuYnZoTmFtZSA9PSBcImhpcFwiKSB7XG4gICAgICAgIG5vZGUuY2hhbm5lbHMgPSBbXCJYcG9zaXRpb25cIiwgXCJZcG9zaXRpb25cIiwgXCJacG9zaXRpb25cIiwgXCJYcm90YXRpb25cIiwgXCJZcm90YXRpb25cIiwgXCJacm90YXRpb25cIl07XG4gICAgICAgIHJldHVybjtcbiAgICB9XG4gICAgbm9kZS5jaGFubmVscyA9IFtdO1xuICAgIGlmICgoX2EgPSBqb2ludCA9PT0gbnVsbCB8fCBqb2ludCA9PT0gdm9pZCAwID8gdm9pZCAwIDogam9pbnQucG9zaXRpb25fa2V5cykgPT09IG51bGwgfHwgX2EgPT09IHZvaWQgMCA/IHZvaWQgMCA6IF9hLmxlbmd0aCkge1xuICAgICAgICBub2RlLmNoYW5uZWxzLnB1c2goXCJYcG9zaXRpb25cIik7XG4gICAgICAgIG5vZGUuY2hhbm5lbHMucHVzaChcIllwb3NpdGlvblwiKTtcbiAgICAgICAgbm9kZS5jaGFubmVscy5wdXNoKFwiWnBvc2l0aW9uXCIpO1xuICAgIH1cbiAgICBub2RlLmNoYW5uZWxzLnB1c2goXCJYcm90YXRpb25cIik7XG4gICAgbm9kZS5jaGFubmVscy5wdXNoKFwiWXJvdGF0aW9uXCIpO1xuICAgIG5vZGUuY2hhbm5lbHMucHVzaChcIlpyb3RhdGlvblwiKTtcbn1cbmZ1bmN0aW9uIGZpbGxLZXlGcmFtZXMoZGF0YSwgYnZoTm9kZSwgZnBzKSB7XG4gICAgY29uc3QgYW5pbUpvaW50cyA9IGRhdGEuam9pbnRzO1xuICAgIGNvbnN0IGxlbmd0aCA9IGV4dHJhY3RGcmFtZXNMZW5ndGgoYW5pbUpvaW50cyk7XG4gICAgY29uc3QgYW5pbVRpbWVzID0gZXh0cmFjdFRpbWVzKGFuaW1Kb2ludHMpO1xuICAgIGNvbnN0IGJ2aFRpbWVzID0gZ2V0VW5pZm9ybVRpbWVzKGRhdGEuZHVyYXRpb24sIDEgLyBmcHMpO1xuICAgIGNvbnN0IGZpeGVkQW5pbVRpbWVzID0gY2xpcFRpbWVzVG9DbG9zZXN0QlZIVGltZShhbmltVGltZXMsIGJ2aFRpbWVzKTtcbiAgICB2aXNpdE5vZGUoYnZoTm9kZSwgKG5vZGUpID0+IHtcbiAgICAgICAgY29uc3Qgam9pbnQgPSBhbmltSm9pbnRzLmZpbmQoKGl0ZW0pID0+IGFsaWFzZXNbaXRlbS5qb2ludF9uYW1lXSA9PT0gbm9kZS5idmhOYW1lKTtcbiAgICAgICAgbm9kZS5vZmZzZXQgPSB7IHg6IDAsIHk6IDAsIHo6IDAgfTtcbiAgICAgICAgaWYgKG5vZGUuYnZoTmFtZSAhPSBcImVuZFwiKSB7XG4gICAgICAgICAgICBmaWxsQ2hhbm5lbHMobm9kZSwgam9pbnQpO1xuICAgICAgICAgICAgbm9kZS5hbmltRnJhbWVzID0gW107XG4gICAgICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IGxlbmd0aDsgaSsrKSB7XG4gICAgICAgICAgICAgICAgbm9kZS5hbmltRnJhbWVzLnB1c2goe1xuICAgICAgICAgICAgICAgICAgICBwb3NpdGlvbjogKGpvaW50ID09PSBudWxsIHx8IGpvaW50ID09PSB2b2lkIDAgPyB2b2lkIDAgOiBqb2ludC5wb3NpdGlvbl9rZXlzW2ldKSB8fCB7IHg6IDAsIHk6IDAsIHo6IDAgfSxcbiAgICAgICAgICAgICAgICAgICAgcm90YXRpb246IChqb2ludCA9PT0gbnVsbCB8fCBqb2ludCA9PT0gdm9pZCAwID8gdm9pZCAwIDogam9pbnQucm90YXRpb25fa2V5c1tpXSkgfHwgeyB4OiAwLCB5OiAwLCB6OiAwIH0sXG4gICAgICAgICAgICAgICAgICAgIHRpbWU6IGFuaW1UaW1lc1tpXVxuICAgICAgICAgICAgICAgIH0pO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgY29uc3QgcG9zaXRpb25zID0gbGVycFZhbHVlcyhub2RlLmFuaW1GcmFtZXMubWFwKChpdGVtKSA9PiBpdGVtLnBvc2l0aW9uKSwgZml4ZWRBbmltVGltZXMsIGJ2aFRpbWVzLCBsZXJwVmVjdG9yKTtcbiAgICAgICAgICAgIGNvbnN0IHJvdGF0aW9ucyA9IGxlcnBWYWx1ZXMobm9kZS5hbmltRnJhbWVzLm1hcCgoaXRlbSkgPT4gdG9RdWF0ZXJuaW9uKGl0ZW0ucm90YXRpb24pKSwgZml4ZWRBbmltVGltZXMsIGJ2aFRpbWVzLCBsZXJwUXVhdGVybmlvbikubWFwKGl0ZW0gPT4gcXVhdGVybmlvblRvRXVsZXJzKGl0ZW0pKTtcbiAgICAgICAgICAgIG5vZGUuYnZoRnJhbWVzID0gW107XG4gICAgICAgICAgICBwb3NpdGlvbnMuZm9yRWFjaCgoaXRlbSwgaSkgPT4gbm9kZS5idmhGcmFtZXMucHVzaCh7XG4gICAgICAgICAgICAgICAgcG9zaXRpb246IGl0ZW0sXG4gICAgICAgICAgICAgICAgcm90YXRpb246IHJvdGF0aW9uc1tpXVxuICAgICAgICAgICAgfSkpO1xuICAgICAgICB9XG4gICAgfSk7XG4gICAgYnZoTm9kZS5idmhUaW1lcyA9IGJ2aFRpbWVzO1xufVxuZnVuY3Rpb24gZ2V0VmFsdWUoYnZoTm9kZSwgY2hhbm5lbCwgZnJhbWVOdW0pIHtcbiAgICBjb25zdCBmcmFtZSA9IGJ2aE5vZGUuYnZoRnJhbWVzW2ZyYW1lTnVtXTtcbiAgICBjb25zdCBrZXkgPSBjaGFubmVsLnRvTG93ZXJDYXNlKClbMF07XG4gICAgY29uc3QgZGF0YSA9IChjaGFubmVsLmluY2x1ZGVzKFwicG9zXCIpKSA/IGZyYW1lLnBvc2l0aW9uIDogZnJhbWUucm90YXRpb247XG4gICAgY29uc3QgdmFsdWUgPSBkYXRhW2tleV07XG4gICAgcmV0dXJuIChNYXRoLmFicyh2YWx1ZSkgPiAwLjAwMDAwMDAxKSA/IHZhbHVlIDogMDtcbn1cbmZ1bmN0aW9uIGdldFZhbHVlcyhidmhOb2RlLCBmcmFtZU51bSkge1xuICAgIHJldHVybiBidmhOb2RlLmNoYW5uZWxzLm1hcCgoaXRlbSkgPT4gZ2V0VmFsdWUoYnZoTm9kZSwgaXRlbSwgZnJhbWVOdW0pKTtcbn1cbmZ1bmN0aW9uIGdldEZyYW1lVmFsdWVzKGJ2aE5vZGUsIGZyYW1lTnVtKSB7XG4gICAgY29uc3QgcmVzdWx0ID0gW107XG4gICAgdmlzaXROb2RlKGJ2aE5vZGUsIChub2RlKSA9PiB7XG4gICAgICAgIGlmICghbm9kZS5jaGFubmVscykge1xuICAgICAgICAgICAgcmV0dXJuO1xuICAgICAgICB9XG4gICAgICAgIHJlc3VsdC51bnNoaWZ0KC4uLmdldFZhbHVlcyhub2RlLCBmcmFtZU51bSkpO1xuICAgIH0sIHRydWUpO1xuICAgIHJldHVybiByZXN1bHQ7XG59XG5mdW5jdGlvbiBnZXRGcmFtZVJvdyhidmhOb2RlLCBmcmFtZU51bSkge1xuICAgIGNvbnN0IHZhbHVlcyA9IGdldEZyYW1lVmFsdWVzKGJ2aE5vZGUsIGZyYW1lTnVtKTtcbiAgICByZXR1cm4gdmFsdWVzLm1hcChpdGVtID0+IGZsb2F0VG9TdHJpbmcoaXRlbSwgNCkpLmpvaW4oXCIgXCIpICsgXCIgXFxuXCI7XG59XG5leHBvcnQgZnVuY3Rpb24gc2VyaWFsaXplQlZIKGJ2aE5vZGUpIHtcbiAgICBsZXQgcmVzdWx0ID0gXCJISUVSQVJDSFlcXG5cIjtcbiAgICByZXN1bHQgKz0gYXBwZW5kTm9kZShidmhOb2RlLCBcIlwiKTtcbiAgICByZXN1bHQgKz0gXCJNT1RJT05cXG5cIjtcbiAgICByZXN1bHQgKz0gXCJGcmFtZXMgXCIgKyBidmhOb2RlLmJ2aFRpbWVzLmxlbmd0aCArIFwiXFxuXCI7XG4gICAgcmVzdWx0ICs9IFwiRnJhbWUgVGltZSBcIiArIGZsb2F0VG9TdHJpbmcoYnZoTm9kZS5idmhUaW1lc1sxXSwgNikgKyBcIlxcblwiO1xuICAgIGZvciAobGV0IGkgPSAwOyBpIDwgYnZoTm9kZS5idmhUaW1lcy5sZW5ndGg7IGkrKykge1xuICAgICAgICByZXN1bHQgKz0gZ2V0RnJhbWVSb3coYnZoTm9kZSwgaSk7XG4gICAgfVxuICAgIHJldHVybiByZXN1bHQ7XG59XG5leHBvcnQgZnVuY3Rpb24gdG9CVkgoZGF0YSkge1xuICAgIGNvbnN0IGJ2aE5vZGUgPSBzdWJUcmVlKGRhdGEuam9pbnRzKTtcbiAgICBmaWxsS2V5RnJhbWVzKGRhdGEsIGJ2aE5vZGUsIDI0KTtcbiAgICByZXR1cm4gYnZoTm9kZTtcbn1cbiIsImV4cG9ydCBjb25zdCBoaWVyYXJjaHkgPSB7XG4gICAgXCJidmhOYW1lXCI6IFwiaGlwXCIsXG4gICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgIHtcbiAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1TcGluZTFcIixcbiAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibVNwaW5lMlwiLFxuICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJhYmRvbWVuXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1TcGluZTNcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibVNwaW5lNFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJjaGVzdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJuZWNrXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImhlYWRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZmlndXJlSGFpclwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRXllUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUV5ZUxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VSb290XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlRXllQWx0UmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlRXllQWx0TGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VGb3JlaGVhZExlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlRm9yZWhlYWRSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VFeWVicm93T3V0ZXJMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZUV5ZWJyb3dDZW50ZXJMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZUV5ZWJyb3dJbm5lckxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlRXllYnJvd091dGVyUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlRXllYnJvd0NlbnRlclJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZUV5ZWJyb3dJbm5lclJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZUV5ZUxpZFVwcGVyTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VFeWVMaWRMb3dlckxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlRXllTGlkVXBwZXJSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VFeWVMaWRMb3dlclJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZUVhcjFMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlRWFyMkxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VFYXIxUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VFYXIyUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VOb3NlTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VOb3NlQ2VudGVyXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZU5vc2VSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VDaGVla0xvd2VyTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VDaGVla1VwcGVyTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VDaGVla0xvd2VyUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlQ2hlZWtVcHBlclJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZUphd1wiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZUNoaW5cIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZVRlZXRoTG93ZXJcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VMaXBMb3dlckxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VMaXBMb3dlclJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlTGlwTG93ZXJDZW50ZXJcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VUb25ndWVCYXNlXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlVG9uZ3VlVGlwXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRmFjZUphd1NoYXBlclwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VGb3JlaGVhZENlbnRlclwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VOb3NlQmFzZVwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VUZWV0aFVwcGVyXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlTGlwVXBwZXJMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VMaXBVcHBlclJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VMaXBDb3JuZXJMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VMaXBDb3JuZXJSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlTGlwVXBwZXJDZW50ZXJcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VFeWVjb3JuZXJJbm5lckxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1GYWNlRXllY29ybmVySW5uZXJSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZhY2VOb3NlQnJpZGdlXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJsQ29sbGFyXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImxTaGxkclwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJsRm9yZUFybVwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJsSGFuZFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZE1pZGRsZTFMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1IYW5kTWlkZGxlMkxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRNaWRkbGUzTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZEluZGV4MUxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRJbmRleDJMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1IYW5kSW5kZXgzTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZFJpbmcxTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZFJpbmcyTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZFJpbmczTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZFBpbmt5MUxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRQaW5reTJMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1IYW5kUGlua3kzTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZFRodW1iMUxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRUaHVtYjJMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1IYW5kVGh1bWIzTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwickNvbGxhclwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJyU2hsZHJcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwickZvcmVBcm1cIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwickhhbmRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRNaWRkbGUxUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRNaWRkbGUyUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRNaWRkbGUzUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRJbmRleDFSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZEluZGV4MlJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1IYW5kSW5kZXgzUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRSaW5nMVJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1IYW5kUmluZzJSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZFJpbmczUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRQaW5reTFSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZFBpbmt5MlJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1IYW5kUGlua3kzUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhhbmRUaHVtYjFSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGFuZFRodW1iMlJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1IYW5kVGh1bWIzUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1XaW5nc1Jvb3RcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibVdpbmcxTGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtV2luZzJMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1XaW5nM0xlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibVdpbmc0TGVmdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1XaW5nNEZhbkxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtV2luZzFSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtV2luZzJSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtV2luZzNSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtV2luZzRSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1XaW5nNEZhblJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICBdXG4gICAgICAgIH0sXG4gICAgICAgIHtcbiAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcInJUaGlnaFwiLFxuICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJyU2hpblwiLFxuICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJyRm9vdFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtRm9vdFJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1Ub2VSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgXVxuICAgICAgICB9LFxuICAgICAgICB7XG4gICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJsVGhpZ2hcIixcbiAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibFNoaW5cIixcbiAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibEZvb3RcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUZvb3RMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1Ub2VMZWZ0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICBdXG4gICAgICAgIH0sXG4gICAgICAgIHtcbiAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1UYWlsMVwiLFxuICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtVGFpbDJcIixcbiAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibVRhaWwzXCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1UYWlsNFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtVGFpbDVcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibVRhaWw2XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcImVuZFwiXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBdXG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICBdXG4gICAgICAgIH0sXG4gICAgICAgIHtcbiAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1Hcm9pblwiLFxuICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJlbmRcIlxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIF1cbiAgICAgICAgfSxcbiAgICAgICAge1xuICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhpbmRMaW1ic1Jvb3RcIixcbiAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhpbmRMaW1iMUxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhpbmRMaW1iMkxlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImNoaWxkcmVuXCI6IFtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhpbmRMaW1iM0xlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiY2hpbGRyZW5cIjogW1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhpbmRMaW1iNExlZnRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGluZExpbWIxUmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwibUhpbmRMaW1iMlJpZ2h0XCIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwiYnZoTmFtZVwiOiBcIm1IaW5kTGltYjNSaWdodFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcImJ2aE5hbWVcIjogXCJtSGluZExpbWI0UmlnaHRcIixcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJjaGlsZHJlblwiOiBbXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCJidmhOYW1lXCI6IFwiZW5kXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIF1cbiAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgXVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIF1cbiAgICAgICAgfVxuICAgIF1cbn07XG4iLCJleHBvcnQgY29uc3QgZmVtYWxlT2Zmc2V0cyA9IHtcbiAgICBcImhpcFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwibUhpbmRMaW1ic1Jvb3RcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDMuMzA3MSxcbiAgICAgICAgXCJ6XCI6IC03Ljg3NFxuICAgIH0sXG4gICAgXCJtSGluZExpbWIxUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTUuMDc4NyxcbiAgICAgICAgXCJ5XCI6IC00LjkyMTMsXG4gICAgICAgIFwielwiOiAtOC4wMzE1XG4gICAgfSxcbiAgICBcIm1IaW5kTGltYjJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjgxMSxcbiAgICAgICAgXCJ5XCI6IC0xOS4zMzA3LFxuICAgICAgICBcInpcIjogMC4wNzg3XG4gICAgfSxcbiAgICBcIm1IaW5kTGltYjNSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjExODEsXG4gICAgICAgIFwieVwiOiAtMTguNDI1MixcbiAgICAgICAgXCJ6XCI6IC0xLjE4MTFcbiAgICB9LFxuICAgIFwibUhpbmRMaW1iNFJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMi40MDE2LFxuICAgICAgICBcInpcIjogNC40MDk0XG4gICAgfSxcbiAgICBcImVuZF9tSGluZExpbWI0UmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuMzE1LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDQuMTMzOVxuICAgIH0sXG4gICAgXCJtSGluZExpbWIxTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiA1LjA3ODcsXG4gICAgICAgIFwieVwiOiAtNC45MjEzLFxuICAgICAgICBcInpcIjogLTguMDMxNVxuICAgIH0sXG4gICAgXCJtSGluZExpbWIyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS44MTEsXG4gICAgICAgIFwieVwiOiAtMTkuMzMwNyxcbiAgICAgICAgXCJ6XCI6IDAuMDc4N1xuICAgIH0sXG4gICAgXCJtSGluZExpbWIzTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC4xMTgxLFxuICAgICAgICBcInlcIjogLTE4LjQyNTIsXG4gICAgICAgIFwielwiOiAtMS4xODExXG4gICAgfSxcbiAgICBcIm1IaW5kTGltYjRMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMi40MDE2LFxuICAgICAgICBcInpcIjogNC40MDk0XG4gICAgfSxcbiAgICBcImVuZF9tSGluZExpbWI0TGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjMxNSxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiA0LjEzMzlcbiAgICB9LFxuICAgIFwibUdyb2luXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMy44MTg5LFxuICAgICAgICBcInpcIjogMi41MTk3XG4gICAgfSxcbiAgICBcImVuZF9tR3JvaW5cIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0yLjU5ODQsXG4gICAgICAgIFwielwiOiAwLjE1NzVcbiAgICB9LFxuICAgIFwibVRhaWwxXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAxLjg1MDQsXG4gICAgICAgIFwielwiOiAtNC41NjY5XG4gICAgfSxcbiAgICBcIm1UYWlsMlwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC03Ljc1NTlcbiAgICB9LFxuICAgIFwibVRhaWwzXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTYuNjE0MlxuICAgIH0sXG4gICAgXCJtVGFpbDRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtNS41OTA1XG4gICAgfSxcbiAgICBcIm1UYWlsNVwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC00LjQwOTRcbiAgICB9LFxuICAgIFwibVRhaWw2XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTMuNzAwOFxuICAgIH0sXG4gICAgXCJlbmRfbVRhaWw2XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTMuNTAzOVxuICAgIH0sXG4gICAgXCJsVGhpZ2hcIjoge1xuICAgICAgICBcInhcIjogNC45OTA3LFxuICAgICAgICBcInlcIjogLTEuNjE0MSxcbiAgICAgICAgXCJ6XCI6IDEuMzI5XG4gICAgfSxcbiAgICBcImxTaGluXCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjc5NCxcbiAgICAgICAgXCJ5XCI6IC0xOS4zMzI4LFxuICAgICAgICBcInpcIjogLTAuMDM0OVxuICAgIH0sXG4gICAgXCJsRm9vdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjA1NDMsXG4gICAgICAgIFwieVwiOiAtMTguNDQyOCxcbiAgICAgICAgXCJ6XCI6IC0xLjEzNzNcbiAgICB9LFxuICAgIFwibUZvb3RMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMi4zODY2LFxuICAgICAgICBcInpcIjogNC40MDc3XG4gICAgfSxcbiAgICBcIm1Ub2VMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogNC4yOTEzXG4gICAgfSxcbiAgICBcImVuZF9tVG9lTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuNzg3NFxuICAgIH0sXG4gICAgXCJyVGhpZ2hcIjoge1xuICAgICAgICBcInhcIjogLTUuMDcxMSxcbiAgICAgICAgXCJ5XCI6IC0xLjYxNzYsXG4gICAgICAgIFwielwiOiAxLjMyMzZcbiAgICB9LFxuICAgIFwiclNoaW5cIjoge1xuICAgICAgICBcInhcIjogMS45MTQ4LFxuICAgICAgICBcInlcIjogLTE5LjMyNzYsXG4gICAgICAgIFwielwiOiAtMC4wMzA3XG4gICAgfSxcbiAgICBcInJGb290XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMTguNDQ0NixcbiAgICAgICAgXCJ6XCI6IC0xLjEzNjZcbiAgICB9LFxuICAgIFwibUZvb3RSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTIuMzg3MyxcbiAgICAgICAgXCJ6XCI6IDQuNDA3N1xuICAgIH0sXG4gICAgXCJtVG9lUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiA0LjI5MTNcbiAgICB9LFxuICAgIFwiZW5kX21Ub2VSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuNzg3NFxuICAgIH0sXG4gICAgXCJtU3BpbmUxXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAzLjMxLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJtU3BpbmUyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMy4zMSxcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwiYWJkb21lblwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMy4zMSxcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwibVNwaW5lM1wiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogOC4wNjYsXG4gICAgICAgIFwielwiOiAtMC42MDVcbiAgICB9LFxuICAgIFwibVNwaW5lNFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTguMDY2LFxuICAgICAgICBcInpcIjogMC42MDVcbiAgICB9LFxuICAgIFwiY2hlc3RcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDguMDY2LFxuICAgICAgICBcInpcIjogLTAuNjA1XG4gICAgfSxcbiAgICBcIm1XaW5nc1Jvb3RcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtMC41NDU3XG4gICAgfSxcbiAgICBcIm1XaW5nMVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0zLjc0MDIsXG4gICAgICAgIFwieVwiOiA3LjEyNixcbiAgICAgICAgXCJ6XCI6IC0zLjg5NzZcbiAgICB9LFxuICAgIFwibVdpbmcyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTYuNjUzNSxcbiAgICAgICAgXCJ5XCI6IDIuNjM3OCxcbiAgICAgICAgXCJ6XCI6IC02LjYxNDJcbiAgICB9LFxuICAgIFwibVdpbmczUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTcuMjA0NyxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtNy4xMjZcbiAgICB9LFxuICAgIFwibVdpbmc0RmFuUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTYuODExLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC02LjczMjNcbiAgICB9LFxuICAgIFwiZW5kX21XaW5nNEZhblJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0yLjQ0MDksXG4gICAgICAgIFwieVwiOiAtNi4yNTk4LFxuICAgICAgICBcInpcIjogLTIuNjc3MlxuICAgIH0sXG4gICAgXCJtV2luZzRSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtNi44MTEsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTYuNzMyM1xuICAgIH0sXG4gICAgXCJlbmRfbVdpbmc0UmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTUuMTk2OCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtNS43NDhcbiAgICB9LFxuICAgIFwibVdpbmcxTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAzLjc0MDIsXG4gICAgICAgIFwieVwiOiA3LjEyNixcbiAgICAgICAgXCJ6XCI6IC0zLjg5NzZcbiAgICB9LFxuICAgIFwibVdpbmcyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiA2LjY1MzUsXG4gICAgICAgIFwieVwiOiAyLjYzNzgsXG4gICAgICAgIFwielwiOiAtNi42MTQyXG4gICAgfSxcbiAgICBcIm1XaW5nM0xlZnRcIjoge1xuICAgICAgICBcInhcIjogNy4yMDQ3LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC03LjEyNlxuICAgIH0sXG4gICAgXCJtV2luZzRGYW5MZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDYuODExLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC02LjczMjNcbiAgICB9LFxuICAgIFwiZW5kX21XaW5nNEZhbkxlZnRcIjoge1xuICAgICAgICBcInhcIjogMi40NDA5LFxuICAgICAgICBcInlcIjogLTYuMjU5OCxcbiAgICAgICAgXCJ6XCI6IC0yLjY3NzJcbiAgICB9LFxuICAgIFwibVdpbmc0TGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiA2LjgxMSxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtNi43MzIzXG4gICAgfSxcbiAgICBcImVuZF9tV2luZzRMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDUuMTk2OCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtNS43NDhcbiAgICB9LFxuICAgIFwickNvbGxhclwiOiB7XG4gICAgICAgIFwieFwiOiAtMi44MzQ2LFxuICAgICAgICBcInlcIjogNi41MTE2LFxuICAgICAgICBcInpcIjogLTAuODE1N1xuICAgIH0sXG4gICAgXCJyU2hsZHJcIjoge1xuICAgICAgICBcInhcIjogLTMuMTI2NyxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcInJGb3JlQXJtXCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xMC45MzU0LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwickhhbmRcIjoge1xuICAgICAgICBcInhcIjogLTkuNTIzNixcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcIm1IYW5kVGh1bWIxUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuMDIzNixcbiAgICAgICAgXCJ5XCI6IDAuMTU3NSxcbiAgICAgICAgXCJ6XCI6IDEuMjIwNVxuICAgIH0sXG4gICAgXCJtSGFuZFRodW1iMlJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjI1OTgsXG4gICAgICAgIFwieVwiOiAtMC4wMzk0LFxuICAgICAgICBcInpcIjogMS4xMDI0XG4gICAgfSxcbiAgICBcIm1IYW5kVGh1bWIzUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuMjIwNSxcbiAgICAgICAgXCJ5XCI6IC0wLjAzOTQsXG4gICAgICAgIFwielwiOiAwLjkwNTVcbiAgICB9LFxuICAgIFwiZW5kX21IYW5kVGh1bWIzUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuOTg0MyxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjU5MDZcbiAgICB9LFxuICAgIFwibUhhbmRQaW5reTFSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMy43NDAxLFxuICAgICAgICBcInlcIjogMC4xMTgxLFxuICAgICAgICBcInpcIjogLTEuMjIwNVxuICAgIH0sXG4gICAgXCJtSGFuZFBpbmt5MlJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjk4NDIsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogLTAuOTQ0OVxuICAgIH0sXG4gICAgXCJtSGFuZFBpbmt5M1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjcwODcsXG4gICAgICAgIFwieVwiOiAtMC4xNTc1LFxuICAgICAgICBcInpcIjogLTAuNTkwNlxuICAgIH0sXG4gICAgXCJlbmRfbUhhbmRQaW5reTNSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC42Mjk5LFxuICAgICAgICBcInlcIjogLTAuMTU3NSxcbiAgICAgICAgXCJ6XCI6IC0wLjUxMThcbiAgICB9LFxuICAgIFwibUhhbmRSaW5nMVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0zLjg5NzYsXG4gICAgICAgIFwieVwiOiAwLjM1NDMsXG4gICAgICAgIFwielwiOiAtMC4zOTM3XG4gICAgfSxcbiAgICBcIm1IYW5kUmluZzJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS40OTYxLFxuICAgICAgICBcInlcIjogLTAuMzE1LFxuICAgICAgICBcInpcIjogLTAuNTExOFxuICAgIH0sXG4gICAgXCJtSGFuZFJpbmczUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuNTc0OCxcbiAgICAgICAgXCJ5XCI6IC0wLjM1NDMsXG4gICAgICAgIFwielwiOiAtMC41MTE4XG4gICAgfSxcbiAgICBcImVuZF9tSGFuZFJpbmczUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuMTAyNCxcbiAgICAgICAgXCJ5XCI6IC0wLjIzNjIsXG4gICAgICAgIFwielwiOiAtMC4zOTM3XG4gICAgfSxcbiAgICBcIm1IYW5kSW5kZXgxUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTMuODE4OSxcbiAgICAgICAgXCJ5XCI6IDAuNTkwNSxcbiAgICAgICAgXCJ6XCI6IDEuNDk2MVxuICAgIH0sXG4gICAgXCJtSGFuZEluZGV4MlJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjQxNzMsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogMC42NjkzXG4gICAgfSxcbiAgICBcIm1IYW5kSW5kZXgzUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuMjU5OCxcbiAgICAgICAgXCJ5XCI6IC0wLjIzNjIsXG4gICAgICAgIFwielwiOiAwLjU1MTJcbiAgICB9LFxuICAgIFwiZW5kX21IYW5kSW5kZXgzUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuOTg0MyxcbiAgICAgICAgXCJ5XCI6IC0wLjE1NzUsXG4gICAgICAgIFwielwiOiAwLjQzMzFcbiAgICB9LFxuICAgIFwibUhhbmRNaWRkbGUxUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTMuOTc2NCxcbiAgICAgICAgXCJ5XCI6IDAuNTkwNSxcbiAgICAgICAgXCJ6XCI6IDAuNTExOFxuICAgIH0sXG4gICAgXCJtSGFuZE1pZGRsZTJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS41NzQ4LFxuICAgICAgICBcInlcIjogLTAuMjM2MixcbiAgICAgICAgXCJ6XCI6IC0wLjAzOTRcbiAgICB9LFxuICAgIFwibUhhbmRNaWRkbGUzUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuOTI5MSxcbiAgICAgICAgXCJ5XCI6IC0wLjMxNSxcbiAgICAgICAgXCJ6XCI6IC0wLjAzOTRcbiAgICB9LFxuICAgIFwiZW5kX21IYW5kTWlkZGxlM1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjI5OTIsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogLTAuMDc4N1xuICAgIH0sXG4gICAgXCJsQ29sbGFyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDIuODIyLFxuICAgICAgICBcInlcIjogNi41MTE2LFxuICAgICAgICBcInpcIjogLTAuODE1N1xuICAgIH0sXG4gICAgXCJsU2hsZHJcIjoge1xuICAgICAgICBcInhcIjogMy4xMTAyLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwibEZvcmVBcm1cIjoge1xuICAgICAgICBcInhcIjogMTAuOTM1NCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcImxIYW5kXCI6IHtcbiAgICAgICAgXCJ4XCI6IDkuNTE2NSxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcIm1IYW5kVGh1bWIxTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjAyMzYsXG4gICAgICAgIFwieVwiOiAwLjE1NzUsXG4gICAgICAgIFwielwiOiAxLjIyMDVcbiAgICB9LFxuICAgIFwibUhhbmRUaHVtYjJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuMjU5OCxcbiAgICAgICAgXCJ5XCI6IC0wLjAzOTQsXG4gICAgICAgIFwielwiOiAxLjEwMjRcbiAgICB9LFxuICAgIFwibUhhbmRUaHVtYjNMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuMjIwNSxcbiAgICAgICAgXCJ5XCI6IC0wLjAzOTQsXG4gICAgICAgIFwielwiOiAwLjkwNTVcbiAgICB9LFxuICAgIFwiZW5kX21IYW5kVGh1bWIzTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjk4NDMsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC41OTA2XG4gICAgfSxcbiAgICBcIm1IYW5kUGlua3kxTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAzLjc0MDIsXG4gICAgICAgIFwieVwiOiAwLjExODEsXG4gICAgICAgIFwielwiOiAtMS4yMjA1XG4gICAgfSxcbiAgICBcIm1IYW5kUGlua3kyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjk4NDIsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogLTAuOTQ0OVxuICAgIH0sXG4gICAgXCJtSGFuZFBpbmt5M0xlZnRcIjoge1xuICAgICAgICBcInhcIjogMC43MDg3LFxuICAgICAgICBcInlcIjogLTAuMTU3NSxcbiAgICAgICAgXCJ6XCI6IC0wLjU5MDZcbiAgICB9LFxuICAgIFwiZW5kX21IYW5kUGlua3kzTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjYyOTksXG4gICAgICAgIFwieVwiOiAtMC4xNTc1LFxuICAgICAgICBcInpcIjogLTAuNTExOFxuICAgIH0sXG4gICAgXCJtSGFuZFJpbmcxTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAzLjg5NzYsXG4gICAgICAgIFwieVwiOiAwLjM1NDMsXG4gICAgICAgIFwielwiOiAtMC4zOTM3XG4gICAgfSxcbiAgICBcIm1IYW5kUmluZzJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuNDk2MSxcbiAgICAgICAgXCJ5XCI6IC0wLjMxNSxcbiAgICAgICAgXCJ6XCI6IC0wLjUxMThcbiAgICB9LFxuICAgIFwibUhhbmRSaW5nM0xlZnRcIjoge1xuICAgICAgICBcInhcIjogMS41NzQ4LFxuICAgICAgICBcInlcIjogLTAuMzU0MyxcbiAgICAgICAgXCJ6XCI6IC0wLjUxMThcbiAgICB9LFxuICAgIFwiZW5kX21IYW5kUmluZzNMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuMTAyNCxcbiAgICAgICAgXCJ5XCI6IC0wLjIzNjIsXG4gICAgICAgIFwielwiOiAtMC4zOTM3XG4gICAgfSxcbiAgICBcIm1IYW5kSW5kZXgxTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAzLjgxODksXG4gICAgICAgIFwieVwiOiAwLjU5MDUsXG4gICAgICAgIFwielwiOiAxLjQ5NjFcbiAgICB9LFxuICAgIFwibUhhbmRJbmRleDJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuNDE3MyxcbiAgICAgICAgXCJ5XCI6IC0wLjIzNjIsXG4gICAgICAgIFwielwiOiAwLjY2OTNcbiAgICB9LFxuICAgIFwibUhhbmRJbmRleDNMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuMjU5OCxcbiAgICAgICAgXCJ5XCI6IC0wLjIzNjIsXG4gICAgICAgIFwielwiOiAwLjU1MTJcbiAgICB9LFxuICAgIFwiZW5kX21IYW5kSW5kZXgzTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjk4NDMsXG4gICAgICAgIFwieVwiOiAtMC4xNTc1LFxuICAgICAgICBcInpcIjogMC40MzMxXG4gICAgfSxcbiAgICBcIm1IYW5kTWlkZGxlMUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMy45NzY0LFxuICAgICAgICBcInlcIjogMC41OTA1LFxuICAgICAgICBcInpcIjogMC41MTE4XG4gICAgfSxcbiAgICBcIm1IYW5kTWlkZGxlMkxlZnRcIjoge1xuICAgICAgICBcInhcIjogMS41NzQ4LFxuICAgICAgICBcInlcIjogLTAuMjM2MixcbiAgICAgICAgXCJ6XCI6IC0wLjAzOTRcbiAgICB9LFxuICAgIFwibUhhbmRNaWRkbGUzTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjkyOTEsXG4gICAgICAgIFwieVwiOiAtMC4zMTUsXG4gICAgICAgIFwielwiOiAtMC4wMzk0XG4gICAgfSxcbiAgICBcImVuZF9tSGFuZE1pZGRsZTNMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuMjk5MixcbiAgICAgICAgXCJ5XCI6IC0wLjIzNjIsXG4gICAgICAgIFwielwiOiAtMC4wNzg3XG4gICAgfSxcbiAgICBcIm5lY2tcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDkuODg2MSxcbiAgICAgICAgXCJ6XCI6IC0wLjM3MDVcbiAgICB9LFxuICAgIFwiaGVhZFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMi45Nzc2LFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJtRmFjZVJvb3RcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDEuNjkzNCxcbiAgICAgICAgXCJ6XCI6IDAuOTEwNFxuICAgIH0sXG4gICAgXCJtRmFjZU5vc2VCcmlkZ2VcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAuNzExOCxcbiAgICAgICAgXCJ6XCI6IDMuMzE0XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZU5vc2VCcmlkZ2VcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAuMzE1LFxuICAgICAgICBcInpcIjogMC41OTA2XG4gICAgfSxcbiAgICBcIm1GYWNlRXllY29ybmVySW5uZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC43MTM1LFxuICAgICAgICBcInlcIjogMS4xMTIyLFxuICAgICAgICBcInpcIjogMi43MzEzXG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUV5ZWNvcm5lcklubmVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjYyOTlcbiAgICB9LFxuICAgIFwibUZhY2VFeWVjb3JuZXJJbm5lckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMC43MTM1LFxuICAgICAgICBcInlcIjogMS4xMTIyLFxuICAgICAgICBcInpcIjogMi43MzEzXG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUV5ZWNvcm5lcklubmVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuNjI5OVxuICAgIH0sXG4gICAgXCJtRmFjZVRlZXRoVXBwZXJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0xLjE0NzgsXG4gICAgICAgIFwielwiOiAwLjcyODNcbiAgICB9LFxuICAgIFwibUZhY2VMaXBVcHBlckNlbnRlclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTAuMDk2NCxcbiAgICAgICAgXCJ6XCI6IDEuNjM4OFxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VMaXBVcHBlckNlbnRlclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMC4wNzg3LFxuICAgICAgICBcInpcIjogMS42OTI5XG4gICAgfSxcbiAgICBcIm1GYWNlTGlwQ29ybmVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMC42OTE5LFxuICAgICAgICBcInlcIjogLTAuMzY0MixcbiAgICAgICAgXCJ6XCI6IDEuMDE5N1xuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VMaXBDb3JuZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMi4wMDc5LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDEuNzcxN1xuICAgIH0sXG4gICAgXCJtRmFjZUxpcENvcm5lckxlZnRcIjoge1xuICAgICAgICBcInhcIjogLTAuNjkxOSxcbiAgICAgICAgXCJ5XCI6IC0wLjM2NDIsXG4gICAgICAgIFwielwiOiAxLjAxOTdcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlTGlwQ29ybmVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAyLjAwNzksXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMS43NzE3XG4gICAgfSxcbiAgICBcIm1GYWNlTGlwVXBwZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTAuMTQ3OCxcbiAgICAgICAgXCJ6XCI6IDEuNjM4OFxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VMaXBVcHBlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjU5MDYsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMS42MTQyXG4gICAgfSxcbiAgICBcIm1GYWNlTGlwVXBwZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMC4xNDc4LFxuICAgICAgICBcInpcIjogMS42Mzg4XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUxpcFVwcGVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjU5MDYsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMS42MTQyXG4gICAgfSxcbiAgICBcIm1GYWNlTm9zZUJhc2VcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0wLjcwNDcsXG4gICAgICAgIFwielwiOiAzLjQyMzJcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlTm9zZUJhc2VcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjU1MTJcbiAgICB9LFxuICAgIFwibUZhY2VGb3JlaGVhZENlbnRlclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMi4zMDI3LFxuICAgICAgICBcInpcIjogMi41MjM3XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUZvcmVoZWFkQ2VudGVyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMS40MTczXG4gICAgfSxcbiAgICBcIm1GYWNlSmF3U2hhcGVyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VKYXdTaGFwZXJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtMC42NjkzXG4gICAgfSxcbiAgICBcIm1GYWNlSmF3XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMC41MzM5LFxuICAgICAgICBcInpcIjogLTAuMDM2NFxuICAgIH0sXG4gICAgXCJtRmFjZVRlZXRoTG93ZXJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0xLjQ4NjksXG4gICAgICAgIFwielwiOiAwLjc2NDhcbiAgICB9LFxuICAgIFwibUZhY2VUb25ndWVCYXNlXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLjE4MjEsXG4gICAgICAgIFwielwiOiAxLjQyMDNcbiAgICB9LFxuICAgIFwibUZhY2VUb25ndWVUaXBcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAuMjU0OSxcbiAgICAgICAgXCJ6XCI6IDAuODAxMlxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VUb25ndWVUaXBcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjM5MzdcbiAgICB9LFxuICAgIFwibUZhY2VMaXBMb3dlckNlbnRlclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDEuNjM4OFxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VMaXBMb3dlckNlbnRlclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMC4wNzg3LFxuICAgICAgICBcInpcIjogMS41NzQ4XG4gICAgfSxcbiAgICBcIm1GYWNlTGlwTG93ZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDEuNjM4OFxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VMaXBMb3dlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjY2OTMsXG4gICAgICAgIFwieVwiOiAwLjE5NjksXG4gICAgICAgIFwielwiOiAxLjMzODZcbiAgICB9LFxuICAgIFwibUZhY2VMaXBMb3dlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAxLjYzODhcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlTGlwTG93ZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuNjY5MyxcbiAgICAgICAgXCJ5XCI6IDAuMTk2OSxcbiAgICAgICAgXCJ6XCI6IDEuMzM4NlxuICAgIH0sXG4gICAgXCJtRmFjZUNoaW5cIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0yLjAyMTIsXG4gICAgICAgIFwielwiOiAyLjUzMVxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VDaGluXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMC43MDg3LFxuICAgICAgICBcInpcIjogMC44MjY4XG4gICAgfSxcbiAgICBcIm1GYWNlQ2hlZWtVcHBlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjI2NjMsXG4gICAgICAgIFwieVwiOiAtMC4yMzEzLFxuICAgICAgICBcInpcIjogMi41NDkyXG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUNoZWVrVXBwZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC41OTA2LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuODY2MVxuICAgIH0sXG4gICAgXCJtRmFjZUNoZWVrTG93ZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS4yNjYzLFxuICAgICAgICBcInlcIjogLTEuMTk0MSxcbiAgICAgICAgXCJ6XCI6IDEuODIwOVxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VDaGVla0xvd2VyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuMTgxMSxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjUxMThcbiAgICB9LFxuICAgIFwibUZhY2VDaGVla1VwcGVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjI2NjMsXG4gICAgICAgIFwieVwiOiAtMC4yMzEzLFxuICAgICAgICBcInpcIjogMi41NDkyXG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUNoZWVrVXBwZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuNTkwNixcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjg2NjFcbiAgICB9LFxuICAgIFwibUZhY2VDaGVla0xvd2VyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjI2NjMsXG4gICAgICAgIFwieVwiOiAtMS4xOTQxLFxuICAgICAgICBcInpcIjogMS44MjA5XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUNoZWVrTG93ZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuMTgxMSxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjUxMThcbiAgICB9LFxuICAgIFwibUZhY2VOb3NlUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuNTU4NyxcbiAgICAgICAgXCJ5XCI6IC0wLjE5NTcsXG4gICAgICAgIFwielwiOiAzLjEzMTlcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlTm9zZVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjE1NzUsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC41OTA2XG4gICAgfSxcbiAgICBcIm1GYWNlTm9zZUNlbnRlclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTAuMDUzNCxcbiAgICAgICAgXCJ6XCI6IDMuNzE0NlxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VOb3NlQ2VudGVyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC45ODQzXG4gICAgfSxcbiAgICBcIm1GYWNlTm9zZUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMC41NTg3LFxuICAgICAgICBcInlcIjogLTAuMTk1NyxcbiAgICAgICAgXCJ6XCI6IDMuMTMxOVxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VOb3NlTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjE1NzUsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC41OTA2XG4gICAgfSxcbiAgICBcIm1GYWNlRWFyMVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0yLjk3OTUsXG4gICAgICAgIFwieVwiOiAwLjA3MTIsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcIm1GYWNlRWFyMlJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjcwODcsXG4gICAgICAgIFwieVwiOiAwLjk4NDIsXG4gICAgICAgIFwielwiOiAtMC43NDhcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRWFyMlJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAxLjI5OTIsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcIm1GYWNlRWFyMUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMi45Nzk1LFxuICAgICAgICBcInlcIjogMC4wNzEyLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJtRmFjZUVhcjJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuNzA4NyxcbiAgICAgICAgXCJ5XCI6IDAuOTg0MixcbiAgICAgICAgXCJ6XCI6IC0wLjc0OFxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFYXIyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMS4yOTkyLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJtRmFjZUV5ZUxpZExvd2VyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuNDE1NixcbiAgICAgICAgXCJ5XCI6IDEuMTgzNCxcbiAgICAgICAgXCJ6XCI6IDIuNjU4NVxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVMaWRMb3dlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMC4yNzU2LFxuICAgICAgICBcInpcIjogMC45NDQ5XG4gICAgfSxcbiAgICBcIm1GYWNlRXllTGlkVXBwZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS40MTU2LFxuICAgICAgICBcInlcIjogMS4xODg3LFxuICAgICAgICBcInpcIjogMi42NTg1XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUV5ZUxpZFVwcGVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAuMTk2OSxcbiAgICAgICAgXCJ6XCI6IDEuMDYzXG4gICAgfSxcbiAgICBcIm1GYWNlRXllTGlkTG93ZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuNDE1NixcbiAgICAgICAgXCJ5XCI6IDEuMTgzNCxcbiAgICAgICAgXCJ6XCI6IDIuNjU4NVxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVMaWRMb3dlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0wLjI3NTYsXG4gICAgICAgIFwielwiOiAwLjk0NDlcbiAgICB9LFxuICAgIFwibUZhY2VFeWVMaWRVcHBlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMS40MTU2LFxuICAgICAgICBcInlcIjogMS4xODg3LFxuICAgICAgICBcInpcIjogMi42NTg1XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUV5ZUxpZFVwcGVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMC4xOTY5LFxuICAgICAgICBcInpcIjogMS4wNjNcbiAgICB9LFxuICAgIFwibUZhY2VFeWVicm93SW5uZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC45MDYsXG4gICAgICAgIFwieVwiOiAxLjg1NzgsXG4gICAgICAgIFwielwiOiAyLjczMTNcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRXllYnJvd0lubmVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAxLjAyMzZcbiAgICB9LFxuICAgIFwibUZhY2VFeWVicm93Q2VudGVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuNjcxMSxcbiAgICAgICAgXCJ5XCI6IDIuMTQ2NyxcbiAgICAgICAgXCJ6XCI6IDIuNTQ5MlxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVicm93Q2VudGVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAxLjA2M1xuICAgIH0sXG4gICAgXCJtRmFjZUV5ZWJyb3dPdXRlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjkxNjgsXG4gICAgICAgIFwieVwiOiAxLjcyOSxcbiAgICAgICAgXCJ6XCI6IDIuMzMwN1xuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVicm93T3V0ZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC41MTE4LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuOTA1NVxuICAgIH0sXG4gICAgXCJtRmFjZUV5ZWJyb3dJbm5lckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMC45MDYsXG4gICAgICAgIFwieVwiOiAxLjg1NzgsXG4gICAgICAgIFwielwiOiAyLjczMTNcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRXllYnJvd0lubmVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDEuMDIzNlxuICAgIH0sXG4gICAgXCJtRmFjZUV5ZWJyb3dDZW50ZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuNjcxMSxcbiAgICAgICAgXCJ5XCI6IDIuMTQ2NyxcbiAgICAgICAgXCJ6XCI6IDIuNTQ5MlxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVicm93Q2VudGVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDEuMDYzXG4gICAgfSxcbiAgICBcIm1GYWNlRXllYnJvd091dGVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjkxNjgsXG4gICAgICAgIFwieVwiOiAxLjcyOSxcbiAgICAgICAgXCJ6XCI6IDIuMzMwN1xuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVicm93T3V0ZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuNTExOCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjkwNTVcbiAgICB9LFxuICAgIFwibUZhY2VGb3JlaGVhZFJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjMwMzUsXG4gICAgICAgIFwieVwiOiAzLjA2MDgsXG4gICAgICAgIFwielwiOiAyLjMzMDdcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRm9yZWhlYWRSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC4xNTc1LFxuICAgICAgICBcInlcIjogMC43MDg3LFxuICAgICAgICBcInpcIjogMC45NDQ5XG4gICAgfSxcbiAgICBcIm1GYWNlRm9yZWhlYWRMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuMzAzNSxcbiAgICAgICAgXCJ5XCI6IDMuMDYwOCxcbiAgICAgICAgXCJ6XCI6IDIuMzMwN1xuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VGb3JlaGVhZExlZnRcIjoge1xuICAgICAgICBcInhcIjogMC4xNTc1LFxuICAgICAgICBcInlcIjogMC43MDg3LFxuICAgICAgICBcInpcIjogMC45NDQ5XG4gICAgfSxcbiAgICBcIm1GYWNlRXllQWx0TGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjQxNTYsXG4gICAgICAgIFwieVwiOiAxLjE4MzYsXG4gICAgICAgIFwielwiOiAyLjY3NTJcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRXllQWx0TGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuOTg0M1xuICAgIH0sXG4gICAgXCJtRmFjZUV5ZUFsdFJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjQxNTYsXG4gICAgICAgIFwieVwiOiAxLjE4MzYsXG4gICAgICAgIFwielwiOiAyLjY3NTRcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRXllQWx0UmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjk4NDNcbiAgICB9LFxuICAgIFwibUV5ZUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMS40MjAzLFxuICAgICAgICBcInlcIjogMi44NzcsXG4gICAgICAgIFwielwiOiAzLjU4NTdcbiAgICB9LFxuICAgIFwiZW5kX21FeWVMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC45ODQzXG4gICAgfSxcbiAgICBcIm1FeWVSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS40MjAzLFxuICAgICAgICBcInlcIjogMi44NzcsXG4gICAgICAgIFwielwiOiAzLjU4NTlcbiAgICB9LFxuICAgIFwiZW5kX21FeWVSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuOTg0M1xuICAgIH0sXG4gICAgXCJmaWd1cmVIYWlyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAyLjYwMzgsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcImVuZF9maWd1cmVIYWlyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAxLjI5OTIsXG4gICAgICAgIFwielwiOiAwXG4gICAgfVxufTtcbmV4cG9ydCBjb25zdCBtYWxlT2Zmc2V0cyA9IHtcbiAgICBcImhpcFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwibUhpbmRMaW1ic1Jvb3RcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDMuMzA3MSxcbiAgICAgICAgXCJ6XCI6IC03Ljg3NFxuICAgIH0sXG4gICAgXCJtSGluZExpbWIxUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTUuMDc4NyxcbiAgICAgICAgXCJ5XCI6IC00LjkyMTMsXG4gICAgICAgIFwielwiOiAtOC4wMzE1XG4gICAgfSxcbiAgICBcIm1IaW5kTGltYjJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjkwMTYsXG4gICAgICAgIFwieVwiOiAtMTkuMzMwNyxcbiAgICAgICAgXCJ6XCI6IDAuMDgyN1xuICAgIH0sXG4gICAgXCJtSGluZExpbWIzUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMC4xMjQsXG4gICAgICAgIFwieVwiOiAtMjAuMjY3NyxcbiAgICAgICAgXCJ6XCI6IC0xLjI0MDJcbiAgICB9LFxuICAgIFwibUhpbmRMaW1iNFJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMi40MDE2LFxuICAgICAgICBcInpcIjogNC40MDk0XG4gICAgfSxcbiAgICBcImVuZF9tSGluZExpbWI0UmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuMzE1LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDQuMTMzOVxuICAgIH0sXG4gICAgXCJtSGluZExpbWIxTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiA1LjA3ODcsXG4gICAgICAgIFwieVwiOiAtNC45MjEzLFxuICAgICAgICBcInpcIjogLTguMDMxNVxuICAgIH0sXG4gICAgXCJtSGluZExpbWIyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS45MDE2LFxuICAgICAgICBcInlcIjogLTE5LjMzMDcsXG4gICAgICAgIFwielwiOiAwLjA4MjdcbiAgICB9LFxuICAgIFwibUhpbmRMaW1iM0xlZnRcIjoge1xuICAgICAgICBcInhcIjogLTAuMTI0LFxuICAgICAgICBcInlcIjogLTIwLjI2NzcsXG4gICAgICAgIFwielwiOiAtMS4yNDAyXG4gICAgfSxcbiAgICBcIm1IaW5kTGltYjRMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMi40MDE2LFxuICAgICAgICBcInpcIjogNC40MDk0XG4gICAgfSxcbiAgICBcImVuZF9tSGluZExpbWI0TGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjMxNSxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiA0LjEzMzlcbiAgICB9LFxuICAgIFwibUdyb2luXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMy44MTg5LFxuICAgICAgICBcInpcIjogMi41MTk3XG4gICAgfSxcbiAgICBcImVuZF9tR3JvaW5cIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0yLjU5ODQsXG4gICAgICAgIFwielwiOiAwLjE1NzVcbiAgICB9LFxuICAgIFwibVRhaWwxXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAxLjg1MDQsXG4gICAgICAgIFwielwiOiAtNC41NjY5XG4gICAgfSxcbiAgICBcIm1UYWlsMlwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC03Ljc1NTlcbiAgICB9LFxuICAgIFwibVRhaWwzXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTYuNjE0MlxuICAgIH0sXG4gICAgXCJtVGFpbDRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtNS41OTA1XG4gICAgfSxcbiAgICBcIm1UYWlsNVwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC00LjQwOTRcbiAgICB9LFxuICAgIFwibVRhaWw2XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTMuNzAwOFxuICAgIH0sXG4gICAgXCJlbmRfbVRhaWw2XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTMuNTAzOVxuICAgIH0sXG4gICAgXCJsVGhpZ2hcIjoge1xuICAgICAgICBcInhcIjogNC45OTA3LFxuICAgICAgICBcInlcIjogLTEuNjE0MSxcbiAgICAgICAgXCJ6XCI6IDEuMzI5XG4gICAgfSxcbiAgICBcImxTaGluXCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjg4MzcsXG4gICAgICAgIFwieVwiOiAtMTkuMzMyOCxcbiAgICAgICAgXCJ6XCI6IC0wLjAzNjdcbiAgICB9LFxuICAgIFwibEZvb3RcIjoge1xuICAgICAgICBcInhcIjogMC4wNTcsXG4gICAgICAgIFwieVwiOiAtMjAuMjg3MSxcbiAgICAgICAgXCJ6XCI6IC0xLjE5NDFcbiAgICB9LFxuICAgIFwibUZvb3RMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMi4zODY2LFxuICAgICAgICBcInpcIjogNC40MDc3XG4gICAgfSxcbiAgICBcIm1Ub2VMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogNC4yOTEzXG4gICAgfSxcbiAgICBcImVuZF9tVG9lTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuNzg3NFxuICAgIH0sXG4gICAgXCJyVGhpZ2hcIjoge1xuICAgICAgICBcInhcIjogLTUuMDcxMSxcbiAgICAgICAgXCJ5XCI6IC0xLjYxNzYsXG4gICAgICAgIFwielwiOiAxLjMyMzZcbiAgICB9LFxuICAgIFwiclNoaW5cIjoge1xuICAgICAgICBcInhcIjogMi4wMTA1LFxuICAgICAgICBcInlcIjogLTE5LjMyNzYsXG4gICAgICAgIFwielwiOiAtMC4wMzIyXG4gICAgfSxcbiAgICBcInJGb290XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMjAuMjg5MSxcbiAgICAgICAgXCJ6XCI6IC0xLjE5MzRcbiAgICB9LFxuICAgIFwibUZvb3RSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTIuMzg3MyxcbiAgICAgICAgXCJ6XCI6IDQuNDA3N1xuICAgIH0sXG4gICAgXCJtVG9lUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiA0LjI5MTNcbiAgICB9LFxuICAgIFwiZW5kX21Ub2VSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuNzg3NFxuICAgIH0sXG4gICAgXCJtU3BpbmUxXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAzLjMxLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJtU3BpbmUyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMy4zMSxcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwiYWJkb21lblwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMy4zMSxcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwibVNwaW5lM1wiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogOC40NjkzLFxuICAgICAgICBcInpcIjogLTAuNjA1XG4gICAgfSxcbiAgICBcIm1TcGluZTRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC04LjQ2OTMsXG4gICAgICAgIFwielwiOiAwLjYwNVxuICAgIH0sXG4gICAgXCJjaGVzdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogOC40NjkzLFxuICAgICAgICBcInpcIjogLTAuNjA1XG4gICAgfSxcbiAgICBcIm1XaW5nc1Jvb3RcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtMC41NzMyXG4gICAgfSxcbiAgICBcIm1XaW5nMVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0zLjc0MDIsXG4gICAgICAgIFwieVwiOiA3LjEyNixcbiAgICAgICAgXCJ6XCI6IC01Ljg2NjFcbiAgICB9LFxuICAgIFwibVdpbmcyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTYuNjUzNSxcbiAgICAgICAgXCJ5XCI6IDIuNjM3OCxcbiAgICAgICAgXCJ6XCI6IC02LjYxNDJcbiAgICB9LFxuICAgIFwibVdpbmczUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTcuMjA0NyxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtNy4xMjZcbiAgICB9LFxuICAgIFwibVdpbmc0RmFuUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTYuODExLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC02LjczMjNcbiAgICB9LFxuICAgIFwiZW5kX21XaW5nNEZhblJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0yLjQ0MDksXG4gICAgICAgIFwieVwiOiAtNi4yNTk4LFxuICAgICAgICBcInpcIjogLTIuNjc3MlxuICAgIH0sXG4gICAgXCJtV2luZzRSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtNi44MTEsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogLTYuNzMyM1xuICAgIH0sXG4gICAgXCJlbmRfbVdpbmc0UmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTUuMTk2OCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtNS43NDhcbiAgICB9LFxuICAgIFwibVdpbmcxTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAzLjc0MDIsXG4gICAgICAgIFwieVwiOiA3LjEyNixcbiAgICAgICAgXCJ6XCI6IC01Ljg2NjFcbiAgICB9LFxuICAgIFwibVdpbmcyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiA2LjY1MzUsXG4gICAgICAgIFwieVwiOiAyLjYzNzgsXG4gICAgICAgIFwielwiOiAtNi42MTQyXG4gICAgfSxcbiAgICBcIm1XaW5nM0xlZnRcIjoge1xuICAgICAgICBcInhcIjogNy4yMDQ3LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC03LjEyNlxuICAgIH0sXG4gICAgXCJtV2luZzRGYW5MZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDYuODExLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IC02LjczMjNcbiAgICB9LFxuICAgIFwiZW5kX21XaW5nNEZhbkxlZnRcIjoge1xuICAgICAgICBcInhcIjogMi40NDA5LFxuICAgICAgICBcInlcIjogLTYuMjU5OCxcbiAgICAgICAgXCJ6XCI6IC0yLjY3NzJcbiAgICB9LFxuICAgIFwibVdpbmc0TGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiA2LjgxMSxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtNi43MzIzXG4gICAgfSxcbiAgICBcImVuZF9tV2luZzRMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDUuMTk2OCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtNS43NDhcbiAgICB9LFxuICAgIFwickNvbGxhclwiOiB7XG4gICAgICAgIFwieFwiOiAtMi45ODIzLFxuICAgICAgICBcInlcIjogNi44MzcyLFxuICAgICAgICBcInpcIjogLTAuODU2OVxuICAgIH0sXG4gICAgXCJyU2hsZHJcIjoge1xuICAgICAgICBcInhcIjogLTQuMzc3NCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcInJGb3JlQXJtXCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xNC4zNTI3LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwickhhbmRcIjoge1xuICAgICAgICBcInhcIjogLTEwLjMzMDcsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJtSGFuZFRodW1iMVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjAyMzYsXG4gICAgICAgIFwieVwiOiAwLjE1NzUsXG4gICAgICAgIFwielwiOiAxLjIyMDVcbiAgICB9LFxuICAgIFwibUhhbmRUaHVtYjJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS4yNTk4LFxuICAgICAgICBcInlcIjogLTAuMDM5NCxcbiAgICAgICAgXCJ6XCI6IDEuMTAyNFxuICAgIH0sXG4gICAgXCJtSGFuZFRodW1iM1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjIyMDUsXG4gICAgICAgIFwieVwiOiAtMC4wMzk0LFxuICAgICAgICBcInpcIjogMC45MDU1XG4gICAgfSxcbiAgICBcImVuZF9tSGFuZFRodW1iM1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjk4NDMsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC41OTA2XG4gICAgfSxcbiAgICBcIm1IYW5kUGlua3kxUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTMuNzQwMixcbiAgICAgICAgXCJ5XCI6IDAuMTE4MSxcbiAgICAgICAgXCJ6XCI6IC0xLjIyMDVcbiAgICB9LFxuICAgIFwibUhhbmRQaW5reTJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC45ODQyLFxuICAgICAgICBcInlcIjogLTAuMjM2MixcbiAgICAgICAgXCJ6XCI6IC0wLjk0NDlcbiAgICB9LFxuICAgIFwibUhhbmRQaW5reTNSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC43MDg3LFxuICAgICAgICBcInlcIjogLTAuMTU3NSxcbiAgICAgICAgXCJ6XCI6IC0wLjU5MDZcbiAgICB9LFxuICAgIFwiZW5kX21IYW5kUGlua3kzUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuNjI5OSxcbiAgICAgICAgXCJ5XCI6IC0wLjE1NzUsXG4gICAgICAgIFwielwiOiAtMC41MTE4XG4gICAgfSxcbiAgICBcIm1IYW5kUmluZzFSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMy44OTc2LFxuICAgICAgICBcInlcIjogMC4zNTQzLFxuICAgICAgICBcInpcIjogLTAuMzkzN1xuICAgIH0sXG4gICAgXCJtSGFuZFJpbmcyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuNDk2MSxcbiAgICAgICAgXCJ5XCI6IC0wLjMxNSxcbiAgICAgICAgXCJ6XCI6IC0wLjUxMThcbiAgICB9LFxuICAgIFwibUhhbmRSaW5nM1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjU3NDgsXG4gICAgICAgIFwieVwiOiAtMC4zNTQzLFxuICAgICAgICBcInpcIjogLTAuNTExOFxuICAgIH0sXG4gICAgXCJlbmRfbUhhbmRSaW5nM1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjEwMjQsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogLTAuMzkzN1xuICAgIH0sXG4gICAgXCJtSGFuZEluZGV4MVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0zLjgxODksXG4gICAgICAgIFwieVwiOiAwLjU5MDUsXG4gICAgICAgIFwielwiOiAxLjQ5NjFcbiAgICB9LFxuICAgIFwibUhhbmRJbmRleDJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS40MTczLFxuICAgICAgICBcInlcIjogLTAuMjM2MixcbiAgICAgICAgXCJ6XCI6IDAuNjY5M1xuICAgIH0sXG4gICAgXCJtSGFuZEluZGV4M1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjI1OTgsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogMC41NTEyXG4gICAgfSxcbiAgICBcImVuZF9tSGFuZEluZGV4M1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjk4NDMsXG4gICAgICAgIFwieVwiOiAtMC4xNTc1LFxuICAgICAgICBcInpcIjogMC40MzMxXG4gICAgfSxcbiAgICBcIm1IYW5kTWlkZGxlMVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0zLjk3NjQsXG4gICAgICAgIFwieVwiOiAwLjU5MDUsXG4gICAgICAgIFwielwiOiAwLjUxMThcbiAgICB9LFxuICAgIFwibUhhbmRNaWRkbGUyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuNTc0OCxcbiAgICAgICAgXCJ5XCI6IC0wLjIzNjIsXG4gICAgICAgIFwielwiOiAtMC4wMzk0XG4gICAgfSxcbiAgICBcIm1IYW5kTWlkZGxlM1JpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjkyOTEsXG4gICAgICAgIFwieVwiOiAtMC4zMTUsXG4gICAgICAgIFwielwiOiAtMC4wMzk0XG4gICAgfSxcbiAgICBcImVuZF9tSGFuZE1pZGRsZTNSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS4yOTkyLFxuICAgICAgICBcInlcIjogLTAuMjM2MixcbiAgICAgICAgXCJ6XCI6IC0wLjA3ODdcbiAgICB9LFxuICAgIFwibENvbGxhclwiOiB7XG4gICAgICAgIFwieFwiOiAyLjk2OSxcbiAgICAgICAgXCJ5XCI6IDYuODM3MixcbiAgICAgICAgXCJ6XCI6IC0wLjg1NjlcbiAgICB9LFxuICAgIFwibFNobGRyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDQuMzU0MyxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcImxGb3JlQXJtXCI6IHtcbiAgICAgICAgXCJ4XCI6IDE0LjM1MjcsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJsSGFuZFwiOiB7XG4gICAgICAgIFwieFwiOiAxMC4zMjI5LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDBcbiAgICB9LFxuICAgIFwibUhhbmRUaHVtYjFMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuMDIzNixcbiAgICAgICAgXCJ5XCI6IDAuMTU3NSxcbiAgICAgICAgXCJ6XCI6IDEuMjIwNVxuICAgIH0sXG4gICAgXCJtSGFuZFRodW1iMkxlZnRcIjoge1xuICAgICAgICBcInhcIjogMS4yNTk4LFxuICAgICAgICBcInlcIjogLTAuMDM5NCxcbiAgICAgICAgXCJ6XCI6IDEuMTAyNFxuICAgIH0sXG4gICAgXCJtSGFuZFRodW1iM0xlZnRcIjoge1xuICAgICAgICBcInhcIjogMS4yMjA1LFxuICAgICAgICBcInlcIjogLTAuMDM5NCxcbiAgICAgICAgXCJ6XCI6IDAuOTA1NVxuICAgIH0sXG4gICAgXCJlbmRfbUhhbmRUaHVtYjNMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuOTg0MyxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjU5MDZcbiAgICB9LFxuICAgIFwibUhhbmRQaW5reTFMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDMuNzQwMixcbiAgICAgICAgXCJ5XCI6IDAuMTE4MSxcbiAgICAgICAgXCJ6XCI6IC0xLjIyMDVcbiAgICB9LFxuICAgIFwibUhhbmRQaW5reTJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuOTg0MixcbiAgICAgICAgXCJ5XCI6IC0wLjIzNjIsXG4gICAgICAgIFwielwiOiAtMC45NDQ5XG4gICAgfSxcbiAgICBcIm1IYW5kUGlua3kzTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjcwODcsXG4gICAgICAgIFwieVwiOiAtMC4xNTc1LFxuICAgICAgICBcInpcIjogLTAuNTkwNlxuICAgIH0sXG4gICAgXCJlbmRfbUhhbmRQaW5reTNMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuNjI5OSxcbiAgICAgICAgXCJ5XCI6IC0wLjE1NzUsXG4gICAgICAgIFwielwiOiAtMC41MTE4XG4gICAgfSxcbiAgICBcIm1IYW5kUmluZzFMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDMuODk3NixcbiAgICAgICAgXCJ5XCI6IDAuMzU0MyxcbiAgICAgICAgXCJ6XCI6IC0wLjM5MzdcbiAgICB9LFxuICAgIFwibUhhbmRSaW5nMkxlZnRcIjoge1xuICAgICAgICBcInhcIjogMS40OTYxLFxuICAgICAgICBcInlcIjogLTAuMzE1LFxuICAgICAgICBcInpcIjogLTAuNTExOFxuICAgIH0sXG4gICAgXCJtSGFuZFJpbmczTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjU3NDgsXG4gICAgICAgIFwieVwiOiAtMC4zNTQzLFxuICAgICAgICBcInpcIjogLTAuNTExOFxuICAgIH0sXG4gICAgXCJlbmRfbUhhbmRSaW5nM0xlZnRcIjoge1xuICAgICAgICBcInhcIjogMS4xMDI0LFxuICAgICAgICBcInlcIjogLTAuMjM2MixcbiAgICAgICAgXCJ6XCI6IC0wLjM5MzdcbiAgICB9LFxuICAgIFwibUhhbmRJbmRleDFMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDMuODE4OSxcbiAgICAgICAgXCJ5XCI6IDAuNTkwNSxcbiAgICAgICAgXCJ6XCI6IDEuNDk2MVxuICAgIH0sXG4gICAgXCJtSGFuZEluZGV4MkxlZnRcIjoge1xuICAgICAgICBcInhcIjogMS40MTczLFxuICAgICAgICBcInlcIjogLTAuMjM2MixcbiAgICAgICAgXCJ6XCI6IDAuNjY5M1xuICAgIH0sXG4gICAgXCJtSGFuZEluZGV4M0xlZnRcIjoge1xuICAgICAgICBcInhcIjogMS4yNTk4LFxuICAgICAgICBcInlcIjogLTAuMjM2MixcbiAgICAgICAgXCJ6XCI6IDAuNTUxMlxuICAgIH0sXG4gICAgXCJlbmRfbUhhbmRJbmRleDNMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuOTg0MyxcbiAgICAgICAgXCJ5XCI6IC0wLjE1NzUsXG4gICAgICAgIFwielwiOiAwLjQzMzFcbiAgICB9LFxuICAgIFwibUhhbmRNaWRkbGUxTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAzLjk3NjQsXG4gICAgICAgIFwieVwiOiAwLjU5MDUsXG4gICAgICAgIFwielwiOiAwLjUxMThcbiAgICB9LFxuICAgIFwibUhhbmRNaWRkbGUyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjU3NDgsXG4gICAgICAgIFwieVwiOiAtMC4yMzYyLFxuICAgICAgICBcInpcIjogLTAuMDM5NFxuICAgIH0sXG4gICAgXCJtSGFuZE1pZGRsZTNMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuOTI5MSxcbiAgICAgICAgXCJ5XCI6IC0wLjMxNSxcbiAgICAgICAgXCJ6XCI6IC0wLjAzOTRcbiAgICB9LFxuICAgIFwiZW5kX21IYW5kTWlkZGxlM0xlZnRcIjoge1xuICAgICAgICBcInhcIjogMS4yOTkyLFxuICAgICAgICBcInlcIjogLTAuMjM2MixcbiAgICAgICAgXCJ6XCI6IC0wLjA3ODdcbiAgICB9LFxuICAgIFwibmVja1wiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMTAuMzgwNCxcbiAgICAgICAgXCJ6XCI6IC0wLjM4OTNcbiAgICB9LFxuICAgIFwiaGVhZFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMy41NzMxLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJtRmFjZVJvb3RcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDEuNjkzNCxcbiAgICAgICAgXCJ6XCI6IDAuOTEwNFxuICAgIH0sXG4gICAgXCJtRmFjZU5vc2VCcmlkZ2VcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAuNzExOCxcbiAgICAgICAgXCJ6XCI6IDMuMzE0XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZU5vc2VCcmlkZ2VcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAuMzE1LFxuICAgICAgICBcInpcIjogMC41OTA2XG4gICAgfSxcbiAgICBcIm1GYWNlRXllY29ybmVySW5uZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC43MTM1LFxuICAgICAgICBcInlcIjogMS4xMTIyLFxuICAgICAgICBcInpcIjogMi43MzEzXG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUV5ZWNvcm5lcklubmVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjYyOTlcbiAgICB9LFxuICAgIFwibUZhY2VFeWVjb3JuZXJJbm5lckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMC43MTM1LFxuICAgICAgICBcInlcIjogMS4xMTIyLFxuICAgICAgICBcInpcIjogMi43MzEzXG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUV5ZWNvcm5lcklubmVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuNjI5OVxuICAgIH0sXG4gICAgXCJtRmFjZVRlZXRoVXBwZXJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0xLjE0NzgsXG4gICAgICAgIFwielwiOiAwLjcyODNcbiAgICB9LFxuICAgIFwibUZhY2VMaXBVcHBlckNlbnRlclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTAuMDk2NCxcbiAgICAgICAgXCJ6XCI6IDEuNjM4OFxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VMaXBVcHBlckNlbnRlclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMC4wNzg3LFxuICAgICAgICBcInpcIjogMS42OTI5XG4gICAgfSxcbiAgICBcIm1GYWNlTGlwQ29ybmVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMC42OTE5LFxuICAgICAgICBcInlcIjogLTAuMzY0MixcbiAgICAgICAgXCJ6XCI6IDEuMDE5N1xuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VMaXBDb3JuZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMi4wMDc5LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDEuNzcxN1xuICAgIH0sXG4gICAgXCJtRmFjZUxpcENvcm5lckxlZnRcIjoge1xuICAgICAgICBcInhcIjogLTAuNjkxOSxcbiAgICAgICAgXCJ5XCI6IC0wLjM2NDIsXG4gICAgICAgIFwielwiOiAxLjAxOTdcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlTGlwQ29ybmVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAyLjAwNzksXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMS43NzE3XG4gICAgfSxcbiAgICBcIm1GYWNlTGlwVXBwZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTAuMTQ3OCxcbiAgICAgICAgXCJ6XCI6IDEuNjM4OFxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VMaXBVcHBlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjU5MDYsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMS42MTQyXG4gICAgfSxcbiAgICBcIm1GYWNlTGlwVXBwZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMC4xNDc4LFxuICAgICAgICBcInpcIjogMS42Mzg4XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUxpcFVwcGVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjU5MDYsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMS42MTQyXG4gICAgfSxcbiAgICBcIm1GYWNlTm9zZUJhc2VcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0wLjcwNDcsXG4gICAgICAgIFwielwiOiAzLjQyMzJcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlTm9zZUJhc2VcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjU1MTJcbiAgICB9LFxuICAgIFwibUZhY2VGb3JlaGVhZENlbnRlclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMi4zMDI3LFxuICAgICAgICBcInpcIjogMi41MjM3XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUZvcmVoZWFkQ2VudGVyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMS40MTczXG4gICAgfSxcbiAgICBcIm1GYWNlSmF3U2hhcGVyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VKYXdTaGFwZXJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAtMC42NjkzXG4gICAgfSxcbiAgICBcIm1GYWNlSmF3XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMC41MzM5LFxuICAgICAgICBcInpcIjogLTAuMDM2NFxuICAgIH0sXG4gICAgXCJtRmFjZVRlZXRoTG93ZXJcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0xLjQ4NjksXG4gICAgICAgIFwielwiOiAwLjc2NDhcbiAgICB9LFxuICAgIFwibUZhY2VUb25ndWVCYXNlXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLjE4MjEsXG4gICAgICAgIFwielwiOiAxLjQyMDNcbiAgICB9LFxuICAgIFwibUZhY2VUb25ndWVUaXBcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAuMjU0OSxcbiAgICAgICAgXCJ6XCI6IDAuODAxMlxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VUb25ndWVUaXBcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjM5MzdcbiAgICB9LFxuICAgIFwibUZhY2VMaXBMb3dlckNlbnRlclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDEuNjM4OFxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VMaXBMb3dlckNlbnRlclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMC4wNzg3LFxuICAgICAgICBcInpcIjogMS41NzQ4XG4gICAgfSxcbiAgICBcIm1GYWNlTGlwTG93ZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDEuNjM4OFxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VMaXBMb3dlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjY2OTMsXG4gICAgICAgIFwieVwiOiAwLjE5NjksXG4gICAgICAgIFwielwiOiAxLjMzODZcbiAgICB9LFxuICAgIFwibUZhY2VMaXBMb3dlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAxLjYzODhcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlTGlwTG93ZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuNjY5MyxcbiAgICAgICAgXCJ5XCI6IDAuMTk2OSxcbiAgICAgICAgXCJ6XCI6IDEuMzM4NlxuICAgIH0sXG4gICAgXCJtRmFjZUNoaW5cIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0yLjAyMTIsXG4gICAgICAgIFwielwiOiAyLjUzMVxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VDaGluXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMC43MDg3LFxuICAgICAgICBcInpcIjogMC44MjY4XG4gICAgfSxcbiAgICBcIm1GYWNlQ2hlZWtVcHBlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjI2NjMsXG4gICAgICAgIFwieVwiOiAtMC4yMzEzLFxuICAgICAgICBcInpcIjogMi41NDkyXG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUNoZWVrVXBwZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC41OTA2LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuODY2MVxuICAgIH0sXG4gICAgXCJtRmFjZUNoZWVrTG93ZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS4yNjYzLFxuICAgICAgICBcInlcIjogLTEuMTk0MSxcbiAgICAgICAgXCJ6XCI6IDEuODIwOVxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VDaGVla0xvd2VyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuMTgxMSxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjUxMThcbiAgICB9LFxuICAgIFwibUZhY2VDaGVla1VwcGVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjI2NjMsXG4gICAgICAgIFwieVwiOiAtMC4yMzEzLFxuICAgICAgICBcInpcIjogMi41NDkyXG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUNoZWVrVXBwZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuNTkwNixcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjg2NjFcbiAgICB9LFxuICAgIFwibUZhY2VDaGVla0xvd2VyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjI2NjMsXG4gICAgICAgIFwieVwiOiAtMS4xOTQxLFxuICAgICAgICBcInpcIjogMS44MjA5XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUNoZWVrTG93ZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuMTgxMSxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjUxMThcbiAgICB9LFxuICAgIFwibUZhY2VOb3NlUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTAuNTU4NyxcbiAgICAgICAgXCJ5XCI6IC0wLjE5NTcsXG4gICAgICAgIFwielwiOiAzLjEzMTlcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlTm9zZVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjE1NzUsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC41OTA2XG4gICAgfSxcbiAgICBcIm1GYWNlTm9zZUNlbnRlclwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogLTAuMDUzNCxcbiAgICAgICAgXCJ6XCI6IDMuNzE0NlxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VOb3NlQ2VudGVyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC45ODQzXG4gICAgfSxcbiAgICBcIm1GYWNlTm9zZUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMC41NTg3LFxuICAgICAgICBcInlcIjogLTAuMTk1NyxcbiAgICAgICAgXCJ6XCI6IDMuMTMxOVxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VOb3NlTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLjE1NzUsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC41OTA2XG4gICAgfSxcbiAgICBcIm1GYWNlRWFyMVJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0yLjk3OTUsXG4gICAgICAgIFwieVwiOiAwLjA3MTIsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcIm1GYWNlRWFyMlJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0wLjcwODcsXG4gICAgICAgIFwieVwiOiAwLjk4NDIsXG4gICAgICAgIFwielwiOiAtMC43NDhcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRWFyMlJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAxLjI5OTIsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcIm1GYWNlRWFyMUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMi45Nzk1LFxuICAgICAgICBcInlcIjogMC4wNzEyLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJtRmFjZUVhcjJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuNzA4NyxcbiAgICAgICAgXCJ5XCI6IDAuOTg0MixcbiAgICAgICAgXCJ6XCI6IC0wLjc0OFxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFYXIyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMS4yOTkyLFxuICAgICAgICBcInpcIjogMFxuICAgIH0sXG4gICAgXCJtRmFjZUV5ZUxpZExvd2VyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuNDE1NixcbiAgICAgICAgXCJ5XCI6IDEuMTgzNCxcbiAgICAgICAgXCJ6XCI6IDIuNjU4NVxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVMaWRMb3dlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAtMC4yNzU2LFxuICAgICAgICBcInpcIjogMC45NDQ5XG4gICAgfSxcbiAgICBcIm1GYWNlRXllTGlkVXBwZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS40MTU2LFxuICAgICAgICBcInlcIjogMS4xODg3LFxuICAgICAgICBcInpcIjogMi42NTg1XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUV5ZUxpZFVwcGVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAuMTk2OSxcbiAgICAgICAgXCJ6XCI6IDEuMDYzXG4gICAgfSxcbiAgICBcIm1GYWNlRXllTGlkTG93ZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuNDE1NixcbiAgICAgICAgXCJ5XCI6IDEuMTgzNCxcbiAgICAgICAgXCJ6XCI6IDIuNjU4NVxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVMaWRMb3dlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IC0wLjI3NTYsXG4gICAgICAgIFwielwiOiAwLjk0NDlcbiAgICB9LFxuICAgIFwibUZhY2VFeWVMaWRVcHBlckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMS40MTU2LFxuICAgICAgICBcInlcIjogMS4xODg3LFxuICAgICAgICBcInpcIjogMi42NTg1XG4gICAgfSxcbiAgICBcImVuZF9tRmFjZUV5ZUxpZFVwcGVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMC4xOTY5LFxuICAgICAgICBcInpcIjogMS4wNjNcbiAgICB9LFxuICAgIFwibUZhY2VFeWVicm93SW5uZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC45MDYsXG4gICAgICAgIFwieVwiOiAxLjg1NzgsXG4gICAgICAgIFwielwiOiAyLjczMTNcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRXllYnJvd0lubmVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAxLjAyMzZcbiAgICB9LFxuICAgIFwibUZhY2VFeWVicm93Q2VudGVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogLTEuNjcxMSxcbiAgICAgICAgXCJ5XCI6IDIuMTQ2NyxcbiAgICAgICAgXCJ6XCI6IDIuNTQ5MlxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVicm93Q2VudGVyUmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAxLjA2M1xuICAgIH0sXG4gICAgXCJtRmFjZUV5ZWJyb3dPdXRlclJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjkxNjgsXG4gICAgICAgIFwieVwiOiAxLjcyOSxcbiAgICAgICAgXCJ6XCI6IDIuMzMwN1xuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVicm93T3V0ZXJSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC41MTE4LFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuOTA1NVxuICAgIH0sXG4gICAgXCJtRmFjZUV5ZWJyb3dJbm5lckxlZnRcIjoge1xuICAgICAgICBcInhcIjogMC45MDYsXG4gICAgICAgIFwieVwiOiAxLjg1NzgsXG4gICAgICAgIFwielwiOiAyLjczMTNcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRXllYnJvd0lubmVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDEuMDIzNlxuICAgIH0sXG4gICAgXCJtRmFjZUV5ZWJyb3dDZW50ZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuNjcxMSxcbiAgICAgICAgXCJ5XCI6IDIuMTQ2NyxcbiAgICAgICAgXCJ6XCI6IDIuNTQ5MlxuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVicm93Q2VudGVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDEuMDYzXG4gICAgfSxcbiAgICBcIm1GYWNlRXllYnJvd091dGVyTGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjkxNjgsXG4gICAgICAgIFwieVwiOiAxLjcyOSxcbiAgICAgICAgXCJ6XCI6IDIuMzMwN1xuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VFeWVicm93T3V0ZXJMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAuNTExOCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjkwNTVcbiAgICB9LFxuICAgIFwibUZhY2VGb3JlaGVhZFJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjMwMzUsXG4gICAgICAgIFwieVwiOiAzLjA2MDgsXG4gICAgICAgIFwielwiOiAyLjMzMDdcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRm9yZWhlYWRSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMC4xNTc1LFxuICAgICAgICBcInlcIjogMC43MDg3LFxuICAgICAgICBcInpcIjogMC45NDQ5XG4gICAgfSxcbiAgICBcIm1GYWNlRm9yZWhlYWRMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDEuMzAzNSxcbiAgICAgICAgXCJ5XCI6IDMuMDYwOCxcbiAgICAgICAgXCJ6XCI6IDIuMzMwN1xuICAgIH0sXG4gICAgXCJlbmRfbUZhY2VGb3JlaGVhZExlZnRcIjoge1xuICAgICAgICBcInhcIjogMC4xNTc1LFxuICAgICAgICBcInlcIjogMC43MDg3LFxuICAgICAgICBcInpcIjogMC45NDQ5XG4gICAgfSxcbiAgICBcIm1GYWNlRXllQWx0TGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAxLjQxNTYsXG4gICAgICAgIFwieVwiOiAxLjE4MzYsXG4gICAgICAgIFwielwiOiAyLjY3NTJcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRXllQWx0TGVmdFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuOTg0M1xuICAgIH0sXG4gICAgXCJtRmFjZUV5ZUFsdFJpZ2h0XCI6IHtcbiAgICAgICAgXCJ4XCI6IC0xLjQxNTYsXG4gICAgICAgIFwieVwiOiAxLjE4MzYsXG4gICAgICAgIFwielwiOiAyLjY3NTRcbiAgICB9LFxuICAgIFwiZW5kX21GYWNlRXllQWx0UmlnaHRcIjoge1xuICAgICAgICBcInhcIjogMCxcbiAgICAgICAgXCJ5XCI6IDAsXG4gICAgICAgIFwielwiOiAwLjk4NDNcbiAgICB9LFxuICAgIFwibUV5ZUxlZnRcIjoge1xuICAgICAgICBcInhcIjogMS40MjAzLFxuICAgICAgICBcInlcIjogMi44NzcsXG4gICAgICAgIFwielwiOiAzLjU4NTdcbiAgICB9LFxuICAgIFwiZW5kX21FeWVMZWZ0XCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAwLFxuICAgICAgICBcInpcIjogMC45ODQzXG4gICAgfSxcbiAgICBcIm1FeWVSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAtMS40MjAzLFxuICAgICAgICBcInlcIjogMi44NzcsXG4gICAgICAgIFwielwiOiAzLjU4NTlcbiAgICB9LFxuICAgIFwiZW5kX21FeWVSaWdodFwiOiB7XG4gICAgICAgIFwieFwiOiAwLFxuICAgICAgICBcInlcIjogMCxcbiAgICAgICAgXCJ6XCI6IDAuOTg0M1xuICAgIH0sXG4gICAgXCJmaWd1cmVIYWlyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAyLjYwMzgsXG4gICAgICAgIFwielwiOiAwXG4gICAgfSxcbiAgICBcImVuZF9maWd1cmVIYWlyXCI6IHtcbiAgICAgICAgXCJ4XCI6IDAsXG4gICAgICAgIFwieVwiOiAxLjI5OTIsXG4gICAgICAgIFwielwiOiAwXG4gICAgfVxufTtcbiIsImltcG9ydCB7IHRvRXVsZXJzLCBkaXN0cmlidXRlVmFsdWUgfSBmcm9tIFwiLi91dGlsc1wiO1xuZXhwb3J0IGZ1bmN0aW9uIHBhcnNlQW5pbShhcnJheUJ1ZmZlcikge1xuICAgIGNvbnN0IHZpZXcgPSBuZXcgRGF0YVZpZXcoYXJyYXlCdWZmZXIpO1xuICAgIGxldCBvZmZzZXQgPSAwO1xuICAgIGZ1bmN0aW9uIHJlYWRVMTYoKSB7XG4gICAgICAgIGNvbnN0IHZhbHVlID0gdmlldy5nZXRVaW50MTYob2Zmc2V0LCB0cnVlKTtcbiAgICAgICAgb2Zmc2V0ICs9IDI7XG4gICAgICAgIHJldHVybiB2YWx1ZTtcbiAgICB9XG4gICAgZnVuY3Rpb24gcmVhZFUzMigpIHtcbiAgICAgICAgY29uc3QgdmFsdWUgPSB2aWV3LmdldFVpbnQxNihvZmZzZXQsIHRydWUpO1xuICAgICAgICBvZmZzZXQgKz0gNDtcbiAgICAgICAgcmV0dXJuIHZhbHVlO1xuICAgIH1cbiAgICBmdW5jdGlvbiByZWFkUzMyKCkge1xuICAgICAgICBjb25zdCB2YWx1ZSA9IHZpZXcuZ2V0SW50MzIob2Zmc2V0LCB0cnVlKTtcbiAgICAgICAgb2Zmc2V0ICs9IDQ7XG4gICAgICAgIHJldHVybiB2YWx1ZTtcbiAgICB9XG4gICAgZnVuY3Rpb24gcmVhZEYzMigpIHtcbiAgICAgICAgY29uc3QgdmFsdWUgPSB2aWV3LmdldEZsb2F0MzIob2Zmc2V0LCB0cnVlKTtcbiAgICAgICAgb2Zmc2V0ICs9IDQ7XG4gICAgICAgIHJldHVybiB2YWx1ZTtcbiAgICB9XG4gICAgZnVuY3Rpb24gcmVhZFN0cmluZygpIHtcbiAgICAgICAgbGV0IHN0ciA9IFwiXCI7XG4gICAgICAgIHdoaWxlIChvZmZzZXQgPCB2aWV3LmJ5dGVMZW5ndGgpIHtcbiAgICAgICAgICAgIGNvbnN0IGNoYXIgPSB2aWV3LmdldFVpbnQ4KG9mZnNldCsrKTtcbiAgICAgICAgICAgIGlmIChjaGFyID09PSAwKVxuICAgICAgICAgICAgICAgIGJyZWFrOyAvLyBOVUxMLXRlcm1pbmF0ZWRcbiAgICAgICAgICAgIHN0ciArPSBTdHJpbmcuZnJvbUNoYXJDb2RlKGNoYXIpO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiBzdHI7XG4gICAgfVxuICAgIGZ1bmN0aW9uIHJlYWRGaXhlZFN0cmluZyhzaXplKSB7XG4gICAgICAgIGxldCBzdHIgPSBcIlwiO1xuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IHNpemU7IGkrKykge1xuICAgICAgICAgICAgY29uc3QgY2hhciA9IHZpZXcuZ2V0VWludDgob2Zmc2V0KyspO1xuICAgICAgICAgICAgaWYgKGNoYXIgIT09IDApXG4gICAgICAgICAgICAgICAgc3RyICs9IFN0cmluZy5mcm9tQ2hhckNvZGUoY2hhcik7XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIHN0cjtcbiAgICB9XG4gICAgY29uc3QgdmVyc2lvbiA9IHJlYWRVMTYoKTtcbiAgICBjb25zdCBzdWJfdmVyc2lvbiA9IHJlYWRVMTYoKTtcbiAgICBjb25zdCBiYXNlX3ByaW9yaXR5ID0gcmVhZFMzMigpO1xuICAgIGNvbnN0IGR1cmF0aW9uID0gcmVhZEYzMigpO1xuICAgIGNvbnN0IGVtb3RlX25hbWUgPSByZWFkU3RyaW5nKCk7XG4gICAgY29uc3QgbG9vcF9pbl9wb2ludCA9IHJlYWRGMzIoKTtcbiAgICBjb25zdCBsb29wX291dF9wb2ludCA9IHJlYWRGMzIoKTtcbiAgICBjb25zdCBsb29wID0gcmVhZFMzMigpO1xuICAgIGNvbnN0IGVhc2VfaW5fZHVyYXRpb24gPSByZWFkRjMyKCk7XG4gICAgY29uc3QgZWFzZV9vdXRfZHVyYXRpb24gPSByZWFkRjMyKCk7XG4gICAgY29uc3QgaGFuZF9wb3NlID0gcmVhZFUzMigpO1xuICAgIGNvbnN0IG51bV9qb2ludHMgPSByZWFkVTMyKCk7XG4gICAgY29uc3Qgam9pbnRzID0gW107XG4gICAgZm9yIChsZXQgaSA9IDA7IGkgPCBudW1fam9pbnRzOyBpKyspIHtcbiAgICAgICAgY29uc3Qgam9pbnRfbmFtZSA9IHJlYWRTdHJpbmcoKTtcbiAgICAgICAgY29uc3Qgam9pbnRfcHJpb3JpdHkgPSByZWFkUzMyKCk7XG4gICAgICAgIGNvbnN0IG51bV9yb3Rfa2V5cyA9IHJlYWRTMzIoKTtcbiAgICAgICAgY29uc3Qgcm90YXRpb25fa2V5cyA9IFtdO1xuICAgICAgICBmb3IgKGxldCBqID0gMDsgaiA8IG51bV9yb3Rfa2V5czsgaisrKSB7XG4gICAgICAgICAgICBjb25zdCB0aW1lID0gcmVhZFUxNigpO1xuICAgICAgICAgICAgY29uc3Qgcm90X3ggPSByZWFkVTE2KCk7XG4gICAgICAgICAgICBjb25zdCByb3RfeSA9IHJlYWRVMTYoKTtcbiAgICAgICAgICAgIGNvbnN0IHJvdF96ID0gcmVhZFUxNigpO1xuICAgICAgICAgICAgcm90YXRpb25fa2V5cy5wdXNoKHsgdGltZTogaW50VG9GbG9hdCh0aW1lLCAwLCBkdXJhdGlvbiksIHg6IGludFRvRmxvYXQocm90X3gsIC0xLCAxKSwgeTogaW50VG9GbG9hdChyb3RfeSwgLTEsIDEpLCB6OiBpbnRUb0Zsb2F0KHJvdF96LCAtMSwgMSkgfSk7XG4gICAgICAgIH1cbiAgICAgICAgY29uc3QgbnVtX3Bvc19rZXlzID0gcmVhZFMzMigpO1xuICAgICAgICBjb25zdCBwb3NpdGlvbl9rZXlzID0gW107XG4gICAgICAgIGZvciAobGV0IGogPSAwOyBqIDwgbnVtX3Bvc19rZXlzOyBqKyspIHtcbiAgICAgICAgICAgIGNvbnN0IHRpbWUgPSByZWFkVTE2KCk7XG4gICAgICAgICAgICBjb25zdCBwb3NfeCA9IHJlYWRVMTYoKTtcbiAgICAgICAgICAgIGNvbnN0IHBvc195ID0gcmVhZFUxNigpO1xuICAgICAgICAgICAgY29uc3QgcG9zX3ogPSByZWFkVTE2KCk7XG4gICAgICAgICAgICBwb3NpdGlvbl9rZXlzLnB1c2goeyB0aW1lOiBpbnRUb0Zsb2F0KHRpbWUsIDAsIGR1cmF0aW9uKSwgeDogaW50VG9GbG9hdChwb3NfeCwgLTUsIDUpLCB5OiBpbnRUb0Zsb2F0KHBvc195LCAtNSwgNSksIHo6IGludFRvRmxvYXQocG9zX3osIC01LCA1KSB9KTtcbiAgICAgICAgfVxuICAgICAgICBqb2ludHMucHVzaCh7IGpvaW50X25hbWUsIGpvaW50X3ByaW9yaXR5LCByb3RhdGlvbl9rZXlzLCBwb3NpdGlvbl9rZXlzIH0pO1xuICAgIH1cbiAgICBjb25zdCBudW1fY29uc3RyYWludHMgPSByZWFkUzMyKCk7XG4gICAgY29uc3QgY29uc3RyYWludHMgPSBbXTtcbiAgICBmb3IgKGxldCBpID0gMDsgaSA8IG51bV9jb25zdHJhaW50czsgaSsrKSB7XG4gICAgICAgIGNvbnN0IGNoYWluX2xlbmd0aCA9IHZpZXcuZ2V0VWludDgob2Zmc2V0KyspO1xuICAgICAgICBjb25zdCBjb25zdHJhaW50X3R5cGUgPSB2aWV3LmdldFVpbnQ4KG9mZnNldCsrKTtcbiAgICAgICAgY29uc3Qgc291cmNlX3ZvbHVtZSA9IHJlYWRGaXhlZFN0cmluZygxNik7XG4gICAgICAgIGNvbnN0IHNvdXJjZV9vZmZzZXQgPSBbcmVhZEYzMigpLCByZWFkRjMyKCksIHJlYWRGMzIoKV07XG4gICAgICAgIGNvbnN0IHRhcmdldF92b2x1bWUgPSByZWFkRml4ZWRTdHJpbmcoMTYpO1xuICAgICAgICBjb25zdCB0YXJnZXRfb2Zmc2V0ID0gW3JlYWRGMzIoKSwgcmVhZEYzMigpLCByZWFkRjMyKCldO1xuICAgICAgICBjb25zdCB0YXJnZXRfZGlyID0gW3JlYWRGMzIoKSwgcmVhZEYzMigpLCByZWFkRjMyKCldO1xuICAgICAgICBjb25zdCBlYXNlX2luX3N0YXJ0ID0gcmVhZEYzMigpO1xuICAgICAgICBjb25zdCBlYXNlX2luX3N0b3AgPSByZWFkRjMyKCk7XG4gICAgICAgIGNvbnN0IGVhc2Vfb3V0X3N0YXJ0ID0gcmVhZEYzMigpO1xuICAgICAgICBjb25zdCBlYXNlX291dF9zdG9wID0gcmVhZEYzMigpO1xuICAgICAgICBjb25zdHJhaW50cy5wdXNoKHtcbiAgICAgICAgICAgIGNoYWluX2xlbmd0aCwgY29uc3RyYWludF90eXBlLCBzb3VyY2Vfdm9sdW1lLCBzb3VyY2Vfb2Zmc2V0LFxuICAgICAgICAgICAgdGFyZ2V0X3ZvbHVtZSwgdGFyZ2V0X29mZnNldCwgdGFyZ2V0X2RpciwgZWFzZV9pbl9zdGFydCwgZWFzZV9pbl9zdG9wLFxuICAgICAgICAgICAgZWFzZV9vdXRfc3RhcnQsIGVhc2Vfb3V0X3N0b3BcbiAgICAgICAgfSk7XG4gICAgfVxuICAgIGpvaW50cy5mb3JFYWNoKChpdGVtKSA9PiBpdGVtLnJvdGF0aW9uX2tleXMuZm9yRWFjaCgocm90KSA9PiB7XG4gICAgICAgIGlmICghaXRlbS5ldWxlcl9rZXlzKSB7XG4gICAgICAgICAgICBpdGVtLmV1bGVyX2tleXMgPSBbXTtcbiAgICAgICAgfVxuICAgICAgICBpdGVtLmV1bGVyX2tleXMucHVzaCh0b0V1bGVycyhyb3QpKTtcbiAgICB9KSk7XG4gICAgcmV0dXJuIHsgdmVyc2lvbiwgc3ViX3ZlcnNpb24sIGR1cmF0aW9uLCBlbW90ZV9uYW1lLCBsb29wLCBqb2ludHMsIGNvbnN0cmFpbnRzIH07XG59XG5mdW5jdGlvbiBpbnRUb0Zsb2F0KHZhbCwgbWluLCBtYXgpIHtcbiAgICBjb25zdCBvbmUgPSAobWF4IC0gbWluKSAvIDY1NTM1LjA7XG4gICAgY29uc3QgcmVzdWx0ID0gbWluICsgdmFsICogb25lO1xuICAgIGlmIChNYXRoLmFicyhyZXN1bHQpIDwgb25lKSB7XG4gICAgICAgIHJldHVybiAwO1xuICAgIH1cbiAgICByZXR1cm4gcmVzdWx0O1xufVxuZnVuY3Rpb24gZW51bWVyYXRlKGNvbnRlbnQsIGtleSwgYWx0ZXIpIHtcbiAgICBsZXQgcmVzdWx0ID0gY29udGVudDtcbiAgICBsZXQgY291bnQgPSAwO1xuICAgIHdoaWxlIChyZXN1bHQuaW5jbHVkZXMoXCJcXFwiXCIgKyBrZXkgKyBcIlxcXCJcIikpIHtcbiAgICAgICAgcmVzdWx0ID0gcmVzdWx0LnJlcGxhY2UoXCJcXFwiXCIgKyBrZXkgKyBcIlxcXCJcIiwgXCJcXFwiXCIgKyBhbHRlciArIGNvdW50ICsgXCJcXFwiXCIpO1xuICAgICAgICBjb3VudCArPSAxO1xuICAgIH1cbiAgICByZXR1cm4gcmVzdWx0O1xufVxuZnVuY3Rpb24gY29tcG9zZShub2RlLCBuYW1lKSB7XG4gICAgY29uc3QgY2hpbGRyZW4gPSBbXTtcbiAgICBjb25zdCByZXN1bHQgPSB7fTtcbiAgICBsZXQgY250cyA9IFtdO1xuICAgIGxldCBqbnRzID0gW107XG4gICAgaWYgKE9iamVjdC5rZXlzKG5vZGUpLmluY2x1ZGVzKFwiRW5kIFNpdGVcIikpIHtcbiAgICAgICAgbm9kZS5qbnQxID0gXCJlbmRcIjtcbiAgICB9XG4gICAgT2JqZWN0LmtleXMobm9kZSkuZm9yRWFjaChpdGVtID0+IHtcbiAgICAgICAgaWYgKGl0ZW0uaW5jbHVkZXMoXCJjbnRcIikpIHtcbiAgICAgICAgICAgIGNudHMucHVzaChpdGVtKTtcbiAgICAgICAgICAgIHJldHVybjtcbiAgICAgICAgfVxuICAgICAgICBpZiAoaXRlbS5pbmNsdWRlcyhcImpudFwiKSkge1xuICAgICAgICAgICAgam50cy5wdXNoKGl0ZW0pO1xuICAgICAgICB9XG4gICAgfSk7XG4gICAgY250cyA9IGNudHMuc29ydCgodmFsMSwgdmFsMikgPT4geyByZXR1cm4gcGFyc2VJbnQodmFsMS5yZXBsYWNlKFwiY250XCIsIFwiXCIpKSAtIHBhcnNlSW50KHZhbDIucmVwbGFjZShcImNudFwiLCBcIlwiKSk7IH0pO1xuICAgIGpudHMgPSBqbnRzLnNvcnQoKHZhbDEsIHZhbDIpID0+IHsgcmV0dXJuIHBhcnNlSW50KHZhbDEucmVwbGFjZShcImpudFwiLCBcIlwiKSkgLSBwYXJzZUludCh2YWwyLnJlcGxhY2UoXCJqbnRcIiwgXCJcIikpOyB9KTtcbiAgICBqbnRzLmZvckVhY2goKGl0ZW0sIGkpID0+IGNoaWxkcmVuLnB1c2goY29tcG9zZShub2RlW2NudHNbaV1dLCBub2RlW2l0ZW1dLnRyaW0oKSkpKTtcbiAgICBpZiAobm9kZS5PRkZTRVQpIHtcbiAgICAgICAgY29uc3Qgb2Zmc2V0ID0gbm9kZS5PRkZTRVQudHJpbSgpLnNwbGl0KFwiIFwiKS5maWx0ZXIoKGl0ZW0pID0+ICFpc05hTihwYXJzZUZsb2F0KGl0ZW0pKSkubWFwKHBhcnNlRmxvYXQpO1xuICAgICAgICByZXN1bHQub2Zmc2V0ID0geyB4OiBvZmZzZXRbMF0sIHk6IG9mZnNldFsxXSwgejogb2Zmc2V0WzJdIH07XG4gICAgfVxuICAgIGlmIChub2RlLkNIQU5ORUxTKSB7XG4gICAgICAgIGNvbnN0IGNoYW5uZWxzID0gbm9kZS5DSEFOTkVMUy50cmltKCkuc3BsaXQoXCIgXCIpLmZpbHRlcigoaXRlbSkgPT4gaXNOYU4ocGFyc2VGbG9hdChpdGVtKSkpO1xuICAgICAgICByZXN1bHQuY2hhbm5lbHMgPSBjaGFubmVscztcbiAgICB9XG4gICAgcmVzdWx0LmJ2aE5hbWUgPSBuYW1lO1xuICAgIGlmIChjaGlsZHJlbi5sZW5ndGgpIHtcbiAgICAgICAgcmVzdWx0LmNoaWxkcmVuID0gY2hpbGRyZW47XG4gICAgfVxuICAgIHJldHVybiByZXN1bHQ7XG59XG5mdW5jdGlvbiBjbGVhbnVwKGRhdGEpIHtcbiAgICBkYXRhLmpudDAgPSBkYXRhLlJPT1Q7XG4gICAgcmV0dXJuIGNvbXBvc2UoZGF0YSwgXCJyb290XCIpO1xufVxuZnVuY3Rpb24gcGFyc2VGcmFtZXMocm93cykge1xuICAgIGNvbnN0IHNwbGl0ZWRSb3dzID0gcm93cy5tYXAoaXRlbSA9PiBpdGVtLnNwbGl0KFwiIFwiKS5tYXAoaXRlbSA9PiBpdGVtLnRyaW0oKSkuZmlsdGVyKGl0ZW0gPT4gISFpdGVtKSk7XG4gICAgcmV0dXJuIHNwbGl0ZWRSb3dzLm1hcChpdGVtID0+IGl0ZW0ubWFwKHBhcnNlRmxvYXQpKTtcbn1cbmZ1bmN0aW9uIHBhcnNlRnJhbWVzUGFydChmcmFtZXNQYXJ0KSB7XG4gICAgY29uc3QgZnJhbWVzUm93cyA9IGZyYW1lc1BhcnQuc3BsaXQoXCJcXG5cIik7XG4gICAgbGV0IHRpbWVJbmRleCA9IC0xO1xuICAgIGZvciAobGV0IGkgPSAwOyBpIDwgZnJhbWVzUm93cy5sZW5ndGg7IGkrKykge1xuICAgICAgICBpZiAoZnJhbWVzUm93c1tpXS50b0xvd2VyQ2FzZSgpLmluY2x1ZGVzKFwidGltZVwiKSkge1xuICAgICAgICAgICAgdGltZUluZGV4ID0gaTtcbiAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICB9XG4gICAgfVxuICAgIGlmICh0aW1lSW5kZXggPCAwKSB7XG4gICAgICAgIHJldHVybiB7IGZyYW1lc0xlbmd0aDogMCwgZnJhbWVEdXJhdGlvbjogMCwgZnJhbWVzOiBbXSB9O1xuICAgIH1cbiAgICBjb25zdCBmcmFtZXNMZW5ndGggPSBwYXJzZUludChmcmFtZXNSb3dzW3RpbWVJbmRleCAtIDFdLnNwbGl0KFwiIFwiKS5tYXAoKGl0ZW0pID0+IGl0ZW0udHJpbSgpKS5maWx0ZXIoKGl0ZW0pID0+ICEhaXRlbSkuZmlsdGVyKChpdGVtKSA9PiAhaXNOYU4oaXRlbSkpWzBdKTtcbiAgICBjb25zdCBmcmFtZUR1cmF0aW9uID0gcGFyc2VGbG9hdChmcmFtZXNSb3dzW3RpbWVJbmRleF0uc3BsaXQoXCIgXCIpLm1hcCgoaXRlbSkgPT4gaXRlbS50cmltKCkpLmZpbHRlcigoaXRlbSkgPT4gISFpdGVtKS5maWx0ZXIoKGl0ZW0pID0+ICFpc05hTihpdGVtKSlbMF0pO1xuICAgIHdoaWxlICghZnJhbWVzUm93c1swXS50b0xvd2VyQ2FzZSgpLmluY2x1ZGVzKFwidGltZVwiKSkge1xuICAgICAgICBmcmFtZXNSb3dzLnNoaWZ0KCk7XG4gICAgfVxuICAgIGZyYW1lc1Jvd3Muc2hpZnQoKTtcbiAgICBjb25zdCBmcmFtZXMgPSBwYXJzZUZyYW1lcyhmcmFtZXNSb3dzKTtcbiAgICByZXR1cm4geyBmcmFtZXNMZW5ndGgsIGZyYW1lRHVyYXRpb24sIGZyYW1lcyB9O1xufVxuZXhwb3J0IGZ1bmN0aW9uIGRpc3RyaWJ1dGVTaW5nbGVGcmFtZShoaWVyYXJjaHksIGZyYW1lKSB7XG4gICAgdmFyIF9hLCBfYjtcbiAgICAoX2EgPSBoaWVyYXJjaHkuY2hpbGRyZW4pID09PSBudWxsIHx8IF9hID09PSB2b2lkIDAgPyB2b2lkIDAgOiBfYS50b1JldmVyc2VkKCkuZm9yRWFjaCgoY2hpbGQpID0+IGRpc3RyaWJ1dGVTaW5nbGVGcmFtZShjaGlsZCwgZnJhbWUpKTtcbiAgICBjb25zdCBwb3NpdGlvbiA9IHsgeDogMCwgeTogMCwgejogMCB9O1xuICAgIGNvbnN0IHJvdGF0aW9uID0geyB4OiAwLCB5OiAwLCB6OiAwIH07XG4gICAgKF9iID0gaGllcmFyY2h5LmNoYW5uZWxzKSA9PT0gbnVsbCB8fCBfYiA9PT0gdm9pZCAwID8gdm9pZCAwIDogX2IudG9SZXZlcnNlZCgpLmZvckVhY2goKGl0ZW0pID0+IHtcbiAgICAgICAgY29uc3QgdmFsdWUgPSBmcmFtZS5wb3AoKTtcbiAgICAgICAgZGlzdHJpYnV0ZVZhbHVlKHBvc2l0aW9uLCByb3RhdGlvbiwgaXRlbSwgdmFsdWUpO1xuICAgIH0pO1xuICAgIGlmICghaGllcmFyY2h5LmJ2aEZyYW1lcykge1xuICAgICAgICBoaWVyYXJjaHkuYnZoRnJhbWVzID0gW107XG4gICAgfVxuICAgIGhpZXJhcmNoeS5idmhGcmFtZXMucHVzaCh7IHBvc2l0aW9uLCByb3RhdGlvbiB9KTtcbn1cbmV4cG9ydCBmdW5jdGlvbiBwYXJzZUJWSCh0ZXh0KSB7XG4gICAgY29uc3QgcGFydHMgPSB0ZXh0LnNwbGl0KFwiTU9USU9OXCIpO1xuICAgIGxldCByZXN1bHQgPSBwYXJ0c1swXS5zcGxpdChcIkhJRVJBUkNIWVwiKVsxXTtcbiAgICBcImFiY2RlZmdoaWprbG1ub3BxcnN0dXZ3eHl6QUJDREVGR0hJSktMTU5PUFFSU1RVVldYWVowMTIzNDU2Nzg5XCIuc3BsaXQoXCJcIikuZm9yRWFjaChpdGVtID0+IHtcbiAgICAgICAgcmVzdWx0ID0gcmVzdWx0LnJlcGxhY2VBbGwoaXRlbSArIFwiXFxuXCIsIGl0ZW0gKyBcIlxcXCIsXFxuXCIpO1xuICAgIH0pO1xuICAgIHJlc3VsdCA9IHJlc3VsdC5yZXBsYWNlQWxsKFwiSk9JTlRcIiwgXCJcXFwiSk9JTlRcXFwiOlxcXCJcIik7XG4gICAgcmVzdWx0ID0gcmVzdWx0LnJlcGxhY2VBbGwoXCJPRkZTRVRcIiwgXCJcXFwiT0ZGU0VUXFxcIjpcXFwiXCIpO1xuICAgIHJlc3VsdCA9IHJlc3VsdC5yZXBsYWNlQWxsKFwiQ0hBTk5FTFNcIiwgXCJcXFwiQ0hBTk5FTFNcXFwiOlxcXCJcIik7XG4gICAgcmVzdWx0ID0gcmVzdWx0LnJlcGxhY2VBbGwoXCJST09UXCIsIFwiXFxcIlJPT1RcXFwiOlxcXCJcIik7XG4gICAgcmVzdWx0ID0gcmVzdWx0LnJlcGxhY2VBbGwoXCJFbmQgU2l0ZVwiLCBcIlxcXCJFbmQgU2l0ZVxcXCI6XFxcIlwiKTtcbiAgICByZXN1bHQgPSByZXN1bHQucmVwbGFjZUFsbChcIntcIiwgXCJcXFwiY29udGVudFxcXCI6IHtcIik7XG4gICAgcmVzdWx0ID0gcmVzdWx0LnNwbGl0KFwifVwiKS5tYXAoaXRlbSA9PiB7XG4gICAgICAgIGlmIChpdGVtLnRyaW0oKS5lbmRzV2l0aChcIixcIikpIHtcbiAgICAgICAgICAgIHJldHVybiBpdGVtLnRyaW0oKSArIFwiXFxcImR1bW15XFxcIjoge31cIjtcbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gaXRlbTtcbiAgICB9KS5qb2luKFwifVwiKTtcbiAgICByZXN1bHQgPSByZXN1bHQuc3BsaXQoXCJ9XCIpLm1hcChpdGVtID0+IHtcbiAgICAgICAgaWYgKGl0ZW0udHJpbSgpLnN0YXJ0c1dpdGgoXCJcXFwiSk9JTlRcXFwiXCIpKSB7XG4gICAgICAgICAgICByZXR1cm4gXCIsXCIgKyBpdGVtLnRyaW0oKTtcbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gaXRlbTtcbiAgICB9KS5qb2luKFwifVwiKTtcbiAgICBsZXQgY291bnQgPSAwO1xuICAgIHJlc3VsdCA9IGVudW1lcmF0ZShyZXN1bHQsIFwiSk9JTlRcIiwgXCJqbnRcIik7XG4gICAgcmVzdWx0ID0gZW51bWVyYXRlKHJlc3VsdCwgXCJjb250ZW50XCIsIFwiY250XCIpO1xuICAgIGNvbnN0IGhpZXJhcmNoeSA9IGNsZWFudXAoSlNPTi5wYXJzZShcIntcIiArIHJlc3VsdCArIFwifVwiKSkuY2hpbGRyZW5bMF07XG4gICAgY29uc3QgYW5pbWF0aW9uID0gcGFyc2VGcmFtZXNQYXJ0KHBhcnRzWzFdKTtcbiAgICBoaWVyYXJjaHkuYnZoVGltZXMgPSBbXTtcbiAgICBmb3IgKGxldCBpID0gMDsgaSA8IGFuaW1hdGlvbi5mcmFtZXNMZW5ndGg7IGkrKykge1xuICAgICAgICBoaWVyYXJjaHkuYnZoVGltZXMucHVzaChhbmltYXRpb24uZnJhbWVEdXJhdGlvbiAqIGkpO1xuICAgIH1cbiAgICBhbmltYXRpb24uZnJhbWVzLmZvckVhY2goaXRlbSA9PiBkaXN0cmlidXRlU2luZ2xlRnJhbWUoaGllcmFyY2h5LCBpdGVtKSk7XG4gICAgcmV0dXJuIGhpZXJhcmNoeTtcbn1cbiIsImltcG9ydCB7IFF1YXRlcm5pb24gfSBmcm9tIFwicXVhdGVybmlvblwiO1xuZXhwb3J0IGNvbnN0IFJBRF9UT19ERUcgPSAxODAuMCAvIE1hdGguUEk7XG5leHBvcnQgZnVuY3Rpb24gdG9RdWF0ZXJuaW9uKHYpIHtcbiAgICBjb25zdCB3U3FyID0gMS4wIC0gKHYueCAqIHYueCArIHYueSAqIHYueSArIHYueiAqIHYueik7XG4gICAgcmV0dXJuIG5ldyBRdWF0ZXJuaW9uKHsgeDogdi54LCB5OiB2LnksIHo6IHYueiwgdzogd1NxciA+IDAgPyBNYXRoLnNxcnQod1NxcikgOiAwIH0pO1xufVxuZXhwb3J0IGZ1bmN0aW9uIHRvRXVsZXJzKHRydW5jYXRlZFF1YW50ZXJpb24pIHtcbiAgICByZXR1cm4gdG9RdWF0ZXJuaW9uKHtcbiAgICAgICAgeDogdHJ1bmNhdGVkUXVhbnRlcmlvbi56LFxuICAgICAgICB5OiB0cnVuY2F0ZWRRdWFudGVyaW9uLngsXG4gICAgICAgIHo6IHRydW5jYXRlZFF1YW50ZXJpb24ueVxuICAgIH0pLnRvRXVsZXIoKS5tYXAoaXRlbTEgPT4gaXRlbTEgKiAxODAgLyBNYXRoLlBJKTtcbn1cbmV4cG9ydCBmdW5jdGlvbiBxdWF0ZXJuaW9uVG9FdWxlcnMocXVhdGVybmlvbikge1xuICAgIGNvbnN0IGV1bGVycyA9IG5ldyBRdWF0ZXJuaW9uKHtcbiAgICAgICAgeDogcXVhdGVybmlvbi56LFxuICAgICAgICB5OiBxdWF0ZXJuaW9uLngsXG4gICAgICAgIHo6IHF1YXRlcm5pb24ueSxcbiAgICAgICAgdzogcXVhdGVybmlvbi53XG4gICAgfSkudG9FdWxlcigpLm1hcChpdGVtMSA9PiBpdGVtMSAqIDE4MCAvIE1hdGguUEkpO1xuICAgIHJldHVybiB7IHg6IGV1bGVyc1swXSwgeTogZXVsZXJzWzFdLCB6OiBldWxlcnNbMl0gfTtcbn1cbmV4cG9ydCBmdW5jdGlvbiBhcHBlbmQodGV4dCwgdGltZXMpIHtcbiAgICBsZXQgcmVzdWx0ID0gXCJcIjtcbiAgICBmb3IgKGxldCBpID0gMDsgaSA8IHRpbWVzOyBpKyspIHtcbiAgICAgICAgcmVzdWx0ICs9IHRleHQ7XG4gICAgfVxuICAgIHJldHVybiByZXN1bHQ7XG59XG5leHBvcnQgZnVuY3Rpb24gZmxvYXRUb1N0cmluZyh2YWx1ZSwgZnJhY3Rpb24pIHtcbiAgICByZXR1cm4gdmFsdWUudG9Mb2NhbGVTdHJpbmcoXCJ1bi1VU1wiLCB7IG1pbmltdW1GcmFjdGlvbkRpZ2l0czogZnJhY3Rpb24gfSkucmVwbGFjZUFsbChcIixcIiwgXCJcIik7XG59XG5mdW5jdGlvbiBsZXJwVmFsdWUoeDEsIHgyLCB0KSB7XG4gICAgaWYgKHQgPiAxKSB7XG4gICAgICAgIHJldHVybiB4MjtcbiAgICB9XG4gICAgaWYgKHQgPCAwKSB7XG4gICAgICAgIHJldHVybiB4MTtcbiAgICB9XG4gICAgcmV0dXJuIHgxICsgKHgyIC0geDEpICogdDtcbn1cbmZ1bmN0aW9uIG9wdGltaXplZEFtaW1hdGlvbkxlbmd0aChkdXJhdGlvbiwgb3JpZ2luYWxGcmFtZUxlbmd0aCkge1xuICAgIGlmIChvcmlnaW5hbEZyYW1lTGVuZ3RoID4gZHVyYXRpb24pIHtcbiAgICAgICAgcmV0dXJuIDI7XG4gICAgfVxuICAgIGNvbnN0IGhvdXJMZW5ndGggPSAzNjAwIC8gb3JpZ2luYWxGcmFtZUxlbmd0aDtcbiAgICBsZXQgYmVzdExlbmd0aCA9IDA7XG4gICAgZm9yIChsZXQgaSA9IDE7IGkgPCBob3VyTGVuZ3RoOyBpKyspIHtcbiAgICAgICAgY29uc3Qgb3B0aW1pemVkRnJhbWVMZW5ndGggPSBkdXJhdGlvbiAvIGk7XG4gICAgICAgIGNvbnN0IGVycm9yID0gTWF0aC5hYnMob3JpZ2luYWxGcmFtZUxlbmd0aCAtIG9wdGltaXplZEZyYW1lTGVuZ3RoKTtcbiAgICAgICAgY29uc3QgcHJldkVycm9yID0gTWF0aC5hYnMob3JpZ2luYWxGcmFtZUxlbmd0aCAtIGJlc3RMZW5ndGgpO1xuICAgICAgICBpZiAoZXJyb3IgPCBwcmV2RXJyb3IpIHtcbiAgICAgICAgICAgIGJlc3RMZW5ndGggPSBvcHRpbWl6ZWRGcmFtZUxlbmd0aDtcbiAgICAgICAgfVxuICAgIH1cbiAgICByZXR1cm4gYmVzdExlbmd0aDtcbn1cbmV4cG9ydCBmdW5jdGlvbiBnZXRVbmlmb3JtVGltZXMoZHVyYXRpb24sIHNpbmdsZUZyYW1lRHVyYXRpb24pIHtcbiAgICBjb25zdCB0aW1lcyA9IFtdO1xuICAgIGNvbnN0IG9wdGltaXplZEZyYW1lRHVyYXRpb24gPSBvcHRpbWl6ZWRBbWltYXRpb25MZW5ndGgoZHVyYXRpb24sIHNpbmdsZUZyYW1lRHVyYXRpb24pO1xuICAgIGNvbnN0IGxlbmd0aCA9IGR1cmF0aW9uIC8gb3B0aW1pemVkRnJhbWVEdXJhdGlvbjtcbiAgICBmb3IgKGxldCBpID0gMDsgaSA8IGxlbmd0aDsgaSsrKSB7XG4gICAgICAgIHRpbWVzLnB1c2goaSAqIG9wdGltaXplZEZyYW1lRHVyYXRpb24pO1xuICAgIH1cbiAgICB0aW1lc1tsZW5ndGggLSAxXSA9IGR1cmF0aW9uO1xuICAgIHJldHVybiB0aW1lcztcbn1cbmZ1bmN0aW9uIGNsb3Nlc3QobGVmdCwgcmlnaHQsIHZhbHVlKSB7XG4gICAgaWYgKE1hdGguYWJzKGxlZnQgLSB2YWx1ZSkgPCBNYXRoLmFicyhyaWdodCAtIHZhbHVlKSkge1xuICAgICAgICByZXR1cm4gbGVmdDtcbiAgICB9XG4gICAgcmV0dXJuIHJpZ2h0O1xufVxuZXhwb3J0IGZ1bmN0aW9uIGNsaXBUaW1lc1RvQ2xvc2VzdEJWSFRpbWUoYW5pbVRpbWVzLCBidmhUaW1lcykge1xuICAgIGNvbnN0IGZpeGVkVGltZXMgPSBhbmltVGltZXMubWFwKChpdGVtKSA9PiBpdGVtKTtcbiAgICBmb3IgKGxldCBpID0gMTsgaSA8IGJ2aFRpbWVzLmxlbmd0aDsgaSsrKSB7XG4gICAgICAgIGZvciAobGV0IGogPSAwOyBqIDwgYW5pbVRpbWVzLmxlbmd0aDsgaisrKSB7XG4gICAgICAgICAgICBjb25zdCBhbmltVGltZSA9IGFuaW1UaW1lc1tqXTtcbiAgICAgICAgICAgIGNvbnN0IGJ2aFRpbWVMZWZ0ID0gYnZoVGltZXNbaSAtIDFdO1xuICAgICAgICAgICAgY29uc3QgYnZoVGltZVJpZ2h0ID0gYnZoVGltZXNbaV07XG4gICAgICAgICAgICBpZiAoKGJ2aFRpbWVMZWZ0IDw9IGFuaW1UaW1lKSAmJiAoYW5pbVRpbWUgPD0gYnZoVGltZVJpZ2h0KSkge1xuICAgICAgICAgICAgICAgIGZpeGVkVGltZXNbal0gPSBjbG9zZXN0KGJ2aFRpbWVMZWZ0LCBidmhUaW1lUmlnaHQsIGFuaW1UaW1lKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgIH1cbiAgICByZXR1cm4gZml4ZWRUaW1lcztcbn1cbmZ1bmN0aW9uIGdldEZhY3RvcnMoYnZoVGltZXMsIGFuaW1UaW1lcykge1xuICAgIHJldHVybiBidmhUaW1lcy5tYXAoKGl0ZW0pID0+IHtcbiAgICAgICAgaWYgKGl0ZW0gPD0gYW5pbVRpbWVzWzBdKSB7XG4gICAgICAgICAgICByZXR1cm4ge1xuICAgICAgICAgICAgICAgIGxlZnRBbmltSW5kZXg6IDAsXG4gICAgICAgICAgICAgICAgcmlnaHRBbmltSW5kZXg6IDEsXG4gICAgICAgICAgICAgICAgZmFjdG9yOiAwXG4gICAgICAgICAgICB9O1xuICAgICAgICB9XG4gICAgICAgIGlmIChpdGVtID49IGFuaW1UaW1lc1thbmltVGltZXMubGVuZ3RoIC0gMV0pIHtcbiAgICAgICAgICAgIHJldHVybiB7XG4gICAgICAgICAgICAgICAgbGVmdEFuaW1JbmRleDogYW5pbVRpbWVzLmxlbmd0aCAtIDIsXG4gICAgICAgICAgICAgICAgcmlnaHRBbmltSW5kZXg6IGFuaW1UaW1lcy5sZW5ndGggLSAxLFxuICAgICAgICAgICAgICAgIGZhY3RvcjogMVxuICAgICAgICAgICAgfTtcbiAgICAgICAgfVxuICAgICAgICBmb3IgKGxldCBpID0gMTsgaSA8IGFuaW1UaW1lcy5sZW5ndGg7IGkrKykge1xuICAgICAgICAgICAgY29uc3QgbGVmdFRpbWUgPSBhbmltVGltZXNbaSAtIDFdO1xuICAgICAgICAgICAgY29uc3QgcmlnaHRUaW1lID0gYW5pbVRpbWVzW2ldO1xuICAgICAgICAgICAgaWYgKChsZWZ0VGltZSA8PSBpdGVtKSAmJiAoaXRlbSA8IHJpZ2h0VGltZSkpIHtcbiAgICAgICAgICAgICAgICBjb25zdCByYW5nZVNpemUgPSByaWdodFRpbWUgLSBsZWZ0VGltZTtcbiAgICAgICAgICAgICAgICByZXR1cm4ge1xuICAgICAgICAgICAgICAgICAgICBsZWZ0QW5pbUluZGV4OiBpIC0gMSxcbiAgICAgICAgICAgICAgICAgICAgcmlnaHRBbmltSW5kZXg6IGksXG4gICAgICAgICAgICAgICAgICAgIGZhY3RvcjogKGl0ZW0gLSBsZWZ0VGltZSkgLyByYW5nZVNpemVcbiAgICAgICAgICAgICAgICB9O1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgfSk7XG59XG5leHBvcnQgZnVuY3Rpb24gZGlzdHJpYnV0ZVZhbHVlKHBvc2l0aW9uLCByb3RhdGlvbiwgY2hhbm5lbCwgdmFsdWUpIHtcbiAgICBjb25zdCBrZXkgPSBjaGFubmVsLnRvTG93ZXJDYXNlKClbMF07XG4gICAgY29uc3QgcmVjaXBpZW50ID0gY2hhbm5lbC5pbmNsdWRlcyhcInBvc1wiKSA/IHBvc2l0aW9uIDogcm90YXRpb247XG4gICAgcmVjaXBpZW50W2tleV0gPSB2YWx1ZTtcbn1cbmV4cG9ydCBmdW5jdGlvbiBsZXJwVmVjdG9yKGxlZnRWYWx1ZSwgcmlnaHRWYWx1ZSwgZmFjdG9yKSB7XG4gICAgcmV0dXJuIHtcbiAgICAgICAgeDogbGVycFZhbHVlKGxlZnRWYWx1ZS54LCByaWdodFZhbHVlLngsIGZhY3RvciksXG4gICAgICAgIHk6IGxlcnBWYWx1ZShsZWZ0VmFsdWUueSwgcmlnaHRWYWx1ZS55LCBmYWN0b3IpLFxuICAgICAgICB6OiBsZXJwVmFsdWUobGVmdFZhbHVlLnosIHJpZ2h0VmFsdWUueiwgZmFjdG9yKSxcbiAgICB9O1xufVxuZXhwb3J0IGZ1bmN0aW9uIGxlcnBRdWF0ZXJuaW9uKGxlZnRWYWx1ZSwgcmlnaHRWYWx1ZSwgZmFjdG9yKSB7XG4gICAgcmV0dXJuIGxlZnRWYWx1ZS5zbGVycChyaWdodFZhbHVlKShmYWN0b3IpO1xufVxuZXhwb3J0IGZ1bmN0aW9uIGxlcnBWYWx1ZXModmFsdWVzLCBhbmltVGltZXMsIHVuaWZvcm1UaW1lcywgbGVycEZ1bmN0aW9uKSB7XG4gICAgcmV0dXJuIGdldEZhY3RvcnModW5pZm9ybVRpbWVzLCBhbmltVGltZXMpLm1hcChpdGVtID0+IHtcbiAgICAgICAgY29uc3QgbGVmdEZyYW1lID0gdmFsdWVzW2l0ZW0ubGVmdEFuaW1JbmRleF07XG4gICAgICAgIGNvbnN0IHJpZ2h0RnJhbWUgPSB2YWx1ZXNbaXRlbS5yaWdodEFuaW1JbmRleF07XG4gICAgICAgIHJldHVybiBsZXJwRnVuY3Rpb24obGVmdEZyYW1lLCByaWdodEZyYW1lLCBpdGVtLmZhY3Rvcik7XG4gICAgfSk7XG59XG4iLCIvLyBUaGUgbW9kdWxlIGNhY2hlXG52YXIgX193ZWJwYWNrX21vZHVsZV9jYWNoZV9fID0ge307XG5cbi8vIFRoZSByZXF1aXJlIGZ1bmN0aW9uXG5mdW5jdGlvbiBfX3dlYnBhY2tfcmVxdWlyZV9fKG1vZHVsZUlkKSB7XG5cdC8vIENoZWNrIGlmIG1vZHVsZSBpcyBpbiBjYWNoZVxuXHR2YXIgY2FjaGVkTW9kdWxlID0gX193ZWJwYWNrX21vZHVsZV9jYWNoZV9fW21vZHVsZUlkXTtcblx0aWYgKGNhY2hlZE1vZHVsZSAhPT0gdW5kZWZpbmVkKSB7XG5cdFx0cmV0dXJuIGNhY2hlZE1vZHVsZS5leHBvcnRzO1xuXHR9XG5cdC8vIENyZWF0ZSBhIG5ldyBtb2R1bGUgKGFuZCBwdXQgaXQgaW50byB0aGUgY2FjaGUpXG5cdHZhciBtb2R1bGUgPSBfX3dlYnBhY2tfbW9kdWxlX2NhY2hlX19bbW9kdWxlSWRdID0ge1xuXHRcdC8vIG5vIG1vZHVsZS5pZCBuZWVkZWRcblx0XHQvLyBubyBtb2R1bGUubG9hZGVkIG5lZWRlZFxuXHRcdGV4cG9ydHM6IHt9XG5cdH07XG5cblx0Ly8gRXhlY3V0ZSB0aGUgbW9kdWxlIGZ1bmN0aW9uXG5cdF9fd2VicGFja19tb2R1bGVzX19bbW9kdWxlSWRdKG1vZHVsZSwgbW9kdWxlLmV4cG9ydHMsIF9fd2VicGFja19yZXF1aXJlX18pO1xuXG5cdC8vIFJldHVybiB0aGUgZXhwb3J0cyBvZiB0aGUgbW9kdWxlXG5cdHJldHVybiBtb2R1bGUuZXhwb3J0cztcbn1cblxuIiwiLy8gZGVmaW5lIGdldHRlciBmdW5jdGlvbnMgZm9yIGhhcm1vbnkgZXhwb3J0c1xuX193ZWJwYWNrX3JlcXVpcmVfXy5kID0gKGV4cG9ydHMsIGRlZmluaXRpb24pID0+IHtcblx0Zm9yKHZhciBrZXkgaW4gZGVmaW5pdGlvbikge1xuXHRcdGlmKF9fd2VicGFja19yZXF1aXJlX18ubyhkZWZpbml0aW9uLCBrZXkpICYmICFfX3dlYnBhY2tfcmVxdWlyZV9fLm8oZXhwb3J0cywga2V5KSkge1xuXHRcdFx0T2JqZWN0LmRlZmluZVByb3BlcnR5KGV4cG9ydHMsIGtleSwgeyBlbnVtZXJhYmxlOiB0cnVlLCBnZXQ6IGRlZmluaXRpb25ba2V5XSB9KTtcblx0XHR9XG5cdH1cbn07IiwiX193ZWJwYWNrX3JlcXVpcmVfXy5vID0gKG9iaiwgcHJvcCkgPT4gKE9iamVjdC5wcm90b3R5cGUuaGFzT3duUHJvcGVydHkuY2FsbChvYmosIHByb3ApKSIsIi8vIGRlZmluZSBfX2VzTW9kdWxlIG9uIGV4cG9ydHNcbl9fd2VicGFja19yZXF1aXJlX18uciA9IChleHBvcnRzKSA9PiB7XG5cdGlmKHR5cGVvZiBTeW1ib2wgIT09ICd1bmRlZmluZWQnICYmIFN5bWJvbC50b1N0cmluZ1RhZykge1xuXHRcdE9iamVjdC5kZWZpbmVQcm9wZXJ0eShleHBvcnRzLCBTeW1ib2wudG9TdHJpbmdUYWcsIHsgdmFsdWU6ICdNb2R1bGUnIH0pO1xuXHR9XG5cdE9iamVjdC5kZWZpbmVQcm9wZXJ0eShleHBvcnRzLCAnX19lc01vZHVsZScsIHsgdmFsdWU6IHRydWUgfSk7XG59OyIsImltcG9ydCB7IHBhcnNlQW5pbSwgcGFyc2VCVkggfSBmcm9tIFwiLi9wYXJzZVwiO1xuaW1wb3J0IHsgdG9CVkgsIHNlcmlhbGl6ZUJWSCwgdmlzaXROb2RlIH0gZnJvbSBcIi4vY29udmVydFwiO1xuaW1wb3J0IHsgbWFsZU9mZnNldHMsIGZlbWFsZU9mZnNldHMgfSBmcm9tIFwiLi9vZmZzZXRzXCI7XG5leHBvcnQgeyBwYXJzZUFuaW0sIHBhcnNlQlZILCB0b0JWSCwgc2VyaWFsaXplQlZILCB2aXNpdE5vZGUsIG1hbGVPZmZzZXRzLCBmZW1hbGVPZmZzZXRzIH07XG4iXSwibmFtZXMiOltdLCJzb3VyY2VSb290IjoiIn0=