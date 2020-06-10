//
//  ZoneList.swift
//  showMe
//
//  Created by Z Chameleon on 27/5/20.
//  Copyright © 2020 Z Cha. All rights reserved.
//

import SwiftUI

struct ZoneList: View {
  @EnvironmentObject var userData: UserData
  
  var body: some View {
    NavigationView {
      List {
        Toggle(isOn: $userData.showSelectedOnly) {
          Text("Selected only")
        }
        
        ForEach(userData.zones) { zone in
          if !self.userData.showSelectedOnly || zone.isSelected {
            NavigationLink(destination: ZoneDetail(zone: zone)) {
              ZoneRow(zone: zone)
            }
          }
        }
      }
      .navigationBarTitle(Text("Zones"))
    }
  }
}

struct ZoneList_Previews: PreviewProvider {
    static var previews: some View {
      ZoneList()
        .environmentObject(UserData())
    }
}
