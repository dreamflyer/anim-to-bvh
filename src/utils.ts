import {Vector3} from "./model";
import {Quaternion} from "quaternion";

export const RAD_TO_DEG: number = 180.0 / Math.PI;

export function toQuaternion(v: Vector3): Quaternion {
	const wSqr: number = 1.0 - (v.x * v.x + v.y * v.y + v.z * v.z);
	
	return new Quaternion({x: v.x, y: v.y, z: v.z, w :wSqr > 0 ? Math.sqrt(wSqr) : 0});
}

export function toEulers(truncatedQuanterion: Vector3) {
	return toQuaternion({
		x: truncatedQuanterion.z,
		y: truncatedQuanterion.x,
		z: truncatedQuanterion.y
	}).toEuler().map(item1 => item1 * 180 / Math.PI);
}

export function quaternionToEulers(quaternion: Quaternion): Vector3 {
	const eulers = new Quaternion({
		x: quaternion.z,
		y: quaternion.x,
		z: quaternion.y,
		w: quaternion.w
	}).toEuler().map(item1 => item1 * 180 / Math.PI);
	
	return {x: eulers[0], y: eulers[1], z: eulers[2]};
}

export function append(text: string, times: number): string {
	let result: string = "";
	
	for(let i = 0; i < times; i++) {
		result += text;
	}
	
	return result;
}

export function floatToString(value: number, fraction: number): string {
	return value.toLocaleString("un-US", {minimumFractionDigits: fraction}).replaceAll(",", "")
}

function lerpValue(x1: number, x2: number, t: number) {
	if(t > 1) {
		return x2;
	}
	
	if(t < 0) {
		return x1;
	}
	
	return x1 + (x2 - x1) * t;
}

function optimizedAmimationLength(duration: number, originalFrameLength: number): number {
	if(originalFrameLength > duration) {
		return 2;
	}
	
	const hourLength: number = 3600 / originalFrameLength;
	
	let bestLength = 0;
	
	for(let i = 1; i < hourLength; i++) {
		const optimizedFrameLength = duration / i;
		
		const error = Math.abs(originalFrameLength - optimizedFrameLength);
		
		const prevError = Math.abs(originalFrameLength - bestLength);
		
		if(error < prevError) {
			bestLength = optimizedFrameLength;
		}
	}
	
	return bestLength;
}

export function getUniformTimes(duration: number, singleFrameDuration: number): number[] {
	const times: number[] = [];
	
	const optimizedFrameDuration = optimizedAmimationLength(duration, singleFrameDuration);
	
	const length = duration / optimizedFrameDuration;
	
	for(let i = 0; i < length; i++) {
		times.push(i * optimizedFrameDuration);
	}
	
	times[length - 1] = duration;
	
	return times;
}

function closest(left: number, right: number, value: number): number {
	if(Math.abs(left - value) < Math.abs(right - value)) {
		return left;
	}
	
	return right;
}

export function clipTimesToClosestBVHTime(animTimes: number[], bvhTimes: number[]): number[] {
	const fixedTimes: number[] = animTimes.map((item: number) => item);
	
	for(let i = 1; i < bvhTimes.length; i++) {
		for(let j = 0; j < animTimes.length; j++) {
			const animTime = animTimes[j];
			
			const bvhTimeLeft = bvhTimes[i - 1];
			const bvhTimeRight = bvhTimes[i];
			
			if((bvhTimeLeft <= animTime) &&  (animTime <= bvhTimeRight)) {
				fixedTimes[j] = closest(bvhTimeLeft, bvhTimeRight, animTime);
			}
		}
	}
	
	return fixedTimes;
}

function getFactors(bvhTimes: number[], animTimes: number[]): any[] {
	return bvhTimes.map((item: number) => {
		if(item <= animTimes[0]) {
			return {
				leftAnimIndex: 0,
				rightAnimIndex: 1,
				factor: 0
			};
		}
		
		if(item >= animTimes[animTimes.length - 1]) {
			return {
				leftAnimIndex: animTimes.length - 2,
				rightAnimIndex: animTimes.length - 1,
				factor: 1
			};
		}
		
		for(let i = 1; i < animTimes.length; i++) {
			const leftTime: number = animTimes[i - 1];
			const rightTime: number = animTimes[i];
			
			if((leftTime <= item) && (item < rightTime)) {
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

export function distributeValue(position: any, rotation: any, channel: string, value: number): void {
	const key: string = channel.toLowerCase()[0];
	
	const recipient = channel.includes("pos") ? position : rotation;
	
	recipient[key] = value;
}

export function lerpVector(leftValue: Vector3, rightValue: Vector3, factor: number): Vector3 {
	return {
		x: lerpValue(leftValue.x, rightValue.x, factor),
		y: lerpValue(leftValue.y, rightValue.y, factor),
		z: lerpValue(leftValue.z, rightValue.z, factor),
	}
}

export function lerpQuaternion(leftValue: Quaternion, rightValue: Quaternion, factor: number): Quaternion {
	return leftValue.slerp(rightValue)(factor);
}

export function lerpValues<T>(values: T[], animTimes: number[], uniformTimes: number[], defaultValue: T, lerpFunction: (leftValue: T, rightValue: T, factor: number) => T): T[] {
	if(!values?.length) {
		return uniformTimes.map(item => defaultValue);
	}
	
	if(values.length == 1) {
		return uniformTimes.map(item => values[0]);
	}
	
	const factors: any[] = getFactors(uniformTimes, animTimes);
	
	return factors.map(item => {
		const leftFrame: T = values[item.leftAnimIndex];
		const rightFrame: T = values[item.rightAnimIndex];
		
		return lerpFunction(leftFrame, rightFrame, item.factor);
	});
}