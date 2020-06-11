//
//  ZoneDetail.swift
//  showMe
//
//  Created by Z Chameleon on 26/5/20.
//  Copyright © 2020 Z Cha. All rights reserved.
//

import SwiftUI

struct ZoneDetail: View {
  @EnvironmentObject var userData: UserData
  var zone: Zone
  
  var zoneIndex: Int {
    userData.zones.firstIndex(where: { $0.id == zone.id })!
  }

  var body: some View {
    VStack(alignment: .leading) {

      VStack {
        CircleImage(image: zone.image)
          .padding(.bottom, -130)
        
        VStack(alignment: .leading) {
          
          HStack {
            Text(zone.name)
              .font(.title)
            
            Button(action: {
              self.userData.zones[self.zoneIndex].isSelected.toggle()
            }) {
              if self.userData.zones[self.zoneIndex].isSelected {
                Image(systemName: "star.fill")
                  .foregroundColor(Color.yellow)
              } else {
                Image(systemName: "star")
                  .foregroundColor(Color.gray)
              }
            }
          }
          .padding(40)
          Divider()
          
          TabView {
            ShowEdges(zone: self.zone)
              .tabItem({ Text("Delaunay Edges") })
              .tag(0)
            ShowVoronoi(zone: self.zone)
              .tabItem({ Text("Voronoi Mesh") })
              .tag(1)
            Text(getJSON(using: zone)!)
              .font(.body)
              .tabItem({ Text("Output") })
              .tag(2)
          }
          
          
        }
      }
      .padding()
      Spacer()
    }
    .navigationBarTitle(Text(zone.name), displayMode: .inline)
    
  }
}

struct ZoneDetail_Previews: PreviewProvider {
    static var previews: some View {
      ZoneDetail(zone: zoneData[0])
        .environmentObject(UserData())
    }
}

// Encode the zone as a JSON string
func getJSON(using zone: Zone) -> String? {
  let encoder = JSONEncoder()
  encoder.outputFormatting = .prettyPrinted
  
  do {
    let data = try encoder.encode(zone)
    return String(data: data, encoding: .utf8)
    
  } catch {
    // Simplest possible error catching
    fatalError("Couldn't create JSON string.")
  }
}
