//
//  ShowEdges.swift
//  showMe
//
//  Created by Z Chameleon on 29/5/20.
//  Copyright © 2020 Z Cha. All rights reserved.
//

import SwiftUI

//let listPoints = zoneData[0].points
//let listEdges = zoneData[0].halfEdges
struct ShowEdges: View {
  @EnvironmentObject var userData: UserData
  var triangulation: Triangulation

  // This draws the triangulation
  var body: some View {
    ZStack {
      GeometryReader { geometry in
        Rectangle()
          .stroke(Color.gray)
        
        Path { path in
          let box = self.triangulation.box

          // Need scaling factors
          let scale: CGFloat = min(geometry.size.width, geometry.size.height) / CGFloat(box.scale)
          
          for (e, _) in self.triangulation.triangles.enumerated() {
            if (e > self.triangulation.halfEdges[e]) {
              let p = self.triangulation.points[self.triangulation.triangles[e]]
              let q = self.triangulation.points[self.triangulation.triangles[nextHalfEdge(edge: e)]] // (edge: e)]]
              
              // Get x & y in point space
              let px = CGFloat(p.x - box.left)
              let py = CGFloat(p.y - box.bottom)
              let qx = CGFloat(q.x - box.left)
              let qy = CGFloat(q.y - box.bottom)
              
              // Scale px & py
              path.move(to: CGPoint(x:scale * px, y:scale * py))
              path.addLine(to: CGPoint(x:scale * qx, y:scale * qy))
            }
          }
        }
        .stroke(Color.purple)
      }
      .padding()
    }
}

struct ShowEdges_Previews: PreviewProvider {
  static var previews: some View {
    ShowEdges(triangulation: triangulationData[0])
      .environmentObject(UserData())  //Text("ZoneShow")
  }
}
}

/*
 
 /* Triangle functions */
 func edgesOf(triangle t: Int) -> Array<Int> { [3 * t, 3 * t + 1, 3 * t + 2] }
 func triangleOf(edge e: Int) -> Int { e / 3 }
 
 /* This needs to create a new array - not just modify an existing one */
 func pointsOf(triangle t:Int, using triangles:Array<Int>) -> Array<Int> {
 return edgesOf(triangle: t).map {triangles[$0]}
 }
 
 //func forEachTriangle(with points:Array<Int>, using delaunay:Delaunator, _ callback:) -> {
 //  for t in delaunay.triangles.length / 3; t++) {
 //    callback(t, pointsOfTriangle(delaunay, t).map(p => points[p]));
 //  }
 //}
 //
 */
