//
//  Draw.swift
//  showMe
//
//  Created by Z Chameleon on 6/6/20.
//  Copyright © 2020 Z Cha. All rights reserved.
//

import Foundation

/* Some Simple helper functions */
func nextHalfEdge(edge e:Int) -> Int { (e % 3 == 2) ? e - 2 : e + 1 }
func prevHalfEdge(edge e:Int) -> Int { (e % 3 == 0) ? e + 2 : e - 1 }


/* Triangle functions */
func edgesOf(triangle t: Int) -> Array<Int> { [3 * t, 3 * t + 1, 3 * t + 2] }
func triangleOf(edge e: Int) -> Int { (e / 3) } // e is an integer
func pointsOf(triangle t:Int, using triangles:Array<Int>) -> Array<Int> {edgesOf(triangle: t).map {triangles[$0]}}

func triangleCentre(triangle t:Int, using points:Array<Point>, using triangles:Array<Int>) -> Point  {
  let vertices:Array<Point> = pointsOf(triangle: t, using: triangles).map({points[$0]})
  return circumCentre(first:  vertices[0],
                      second: vertices[1],
                      third:  vertices[2])
}
