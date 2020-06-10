//
//  ZoneRow.swift
//  showMe
//
//  Created by Z Chameleon on 26/5/20.
//  Copyright © 2020 Z Cha. All rights reserved.
//

import SwiftUI

struct ZoneRow: View {
  var zone: Zone
    var body: some View {
      HStack {
        zone.image
          .resizable()
          .frame(width: 50, height: 50)
        Text(zone.name)
        Spacer()
        
        // User selection
        if zone.isSelected {
          Image(systemName: "star.fill")
            .imageScale(.medium)
            .foregroundColor(.yellow)
        }
      }
    }
}

struct ZoneRow_Previews: PreviewProvider {
    static var previews: some View {
      Group {
          ZoneRow(zone: zoneData[0])
          ZoneRow(zone: zoneData[1])
          ZoneRow(zone: zoneData[2])
      }
      .previewLayout(.fixed(width: 300, height: 70))
  }
}
