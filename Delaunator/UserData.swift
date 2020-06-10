//
//  UserData.swift
//  showMe
//
//  Created by Z Chameleon on 28/5/20.
//  Copyright © 2020 Z Cha. All rights reserved.
//

import SwiftUI
import Combine


final class UserData: ObservableObject  {
  @Published var showSelectedOnly = false
  @Published var zones = zoneData
  @Published var triangulations = triangulationData

}
