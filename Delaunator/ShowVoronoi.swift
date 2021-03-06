//
//  ShowVoronoi.swift
//  Delaunator
//
//  Created by Z Chameleon on 29/5/20.
//  Copyright © 2020 Z Cha. All rights reserved.
//

import SwiftUI

struct ShowVoronoi: View {
  @EnvironmentObject var userData: UserData
  var zone: Zone

  // This draws the triangulation
  var body: some View {
    ZStack {
      GeometryReader { geometry in
        
        Rectangle()
          .stroke(Color.gray)
        
        Path { path in
          let triangulation = self.zone.triangulation ?? Triangulation()
          let box = self.zone.box ?? Box()
          
          // Need scaling factors
          let scale: CGFloat = min(geometry.size.width, geometry.size.height) / CGFloat(box.scale)
          
          for (e, _) in triangulation.triangles.enumerated() {
            
            if (e < triangulation.halfEdges[e]) {
              let p = triangleCentre(triangle: triangleOf(edge: e),
                                     using: triangulation.points,
                                     using: triangulation.triangles)
              let q = triangleCentre(triangle: triangleOf(edge: triangulation.halfEdges[e]),
                                     using: triangulation.points,
                                     using: triangulation.triangles)
              
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
        .stroke(Color.blue)
        .clipped(antialiased: false)
      }
    }
    .padding()
  
  }
  
  struct ShowVoronoi_Previews: PreviewProvider {
    static var previews: some View {
      ShowVoronoi(zone: zoneData[0])
        .environmentObject(UserData())  //Text("ZoneShow")
    }
  }
}


